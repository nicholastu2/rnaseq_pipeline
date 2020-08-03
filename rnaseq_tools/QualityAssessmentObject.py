import os
import subprocess
import re
import configparser
import pandas as pd
from glob import glob
import sys
from rnaseq_tools import utils
from rnaseq_tools.StandardDataObject import StandardData
from rnaseq_tools.DatabaseObject import DatabaseObject
from rnaseq_tools.OrganismDataObject import OrganismData
import abc

# turn off SettingWithCopyWarning in pandas
pd.options.mode.chained_assignment = None


class QualityAssessmentObject(StandardData):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = ['bam_file_list', 'count_file_list', 'novoalign_log_list',
                                         'coverage_check_flag','query_path', 'standardized_database_df', 'qual_assess_dir_path']
        # This is a method of adding expected attributes to StandardData from StandardData children
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(QualityAssessmentObject, self).__init__(self._add_expected_attributes, **kwargs)
        # set standardDirectory structure
        self.standardDirectoryStructure()
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'QualityAssessmentObject'

        # if query_path passed, read in as dataframe
        try:
            self.query_df = utils.readInDataframe(self.query_path)
            try:
                if not len(self.query_df[self.query_df.fastqFileName.isnull()]) == 0:
                    raise ValueError('EmptyFastqFileNamesIn %s' %self.query_path)
            except ValueError:
                self.logger.critical('Removing rows with null fastqFileNames from %s' %self.query_path)
                self.query_df = self.query_df[~self.query_df.fastqFileName.isnull()]
        except ValueError:
            self.logger.critical('query_path not valid')
        except FileNotFoundError:
            self.logger.critical('%s  --> query_path not valid' % self.query_path)
        except AttributeError:
            pass

    def compileAlignCountMetadata(self):  # TODO: clean up this, parseAlignmentLog and parseCountFile
        """
        get a list of the filenames in the run_#### file that correspond to a given type
        :returns: a dataframe containing the files according to their suffix
        """
        # instantiate dataframe
        align_df = pd.DataFrame()
        htseq_count_df = pd.DataFrame()
        # extract metadata from novoalign log files
        for log_file in self.novoalign_log_list:
            # extract fastq filename
            fastq_basename = utils.pathBaseName(log_file).replace('_novoalign', '')
            # set sample name in library_metadata_dict
            library_metadata_dict = {"FASTQFILENAME": fastq_basename}
            print('...extracting information from novoalign log for %s' % fastq_basename)
            library_metadata_dict.update(self.parseAlignmentLog(log_file))
            align_df = align_df.append(pd.Series(library_metadata_dict), ignore_index=True)
        print('\nDone parsing novoalign logs\n')
        # extract metadata from count files
        for count_file in self.count_file_list:
            # extract fastq filename
            fastq_basename = utils.pathBaseName(count_file).replace('_read_count', '')
            # set sample name in library_metadata_dict
            library_metadata_dict = {"FASTQFILENAME": fastq_basename}
            print('...extracting count information from count file for %s' % fastq_basename)
            library_metadata_dict.update(self.parseGeneCount(count_file))
            library_metadata_dict['AMBIGUOUS_UNIQUE_PROTEIN_CODING_READS'] = self.uniqueAmbiguousProteinCodingCount(
                fastq_basename)
            htseq_count_df = htseq_count_df.append(pd.Series(library_metadata_dict), ignore_index=True)
        print('\nDone parsing count files\n')
        # concat df_list dataframes together on the common column. CREDIT: https://stackoverflow.com/a/56324303/9708266
        qual_assess_df = pd.merge(align_df, htseq_count_df, on='FASTQFILENAME')
        print('Quantifying noncoding rRNA (rRNA, tRNA and ncRNA)')
        # extract rRNA, tRNA and ncRNA quantification for crypto from bam files -- this takes a long time
        ncRNA_df = self.quantifyNonCodingRna(qual_assess_df)
        # merge this into the qual_assess_df
        qual_assess_df = pd.merge(qual_assess_df, ncRNA_df, on='FASTQFILENAME')
        print('Quantifying intergenic coverage')
        qual_assess_df = self.calculateIntergenicCoverage(qual_assess_df)
        # if coverage_check_flag true, check coverage of perturbed genes
        try:
            if self.coverage_check_flag:
                coverage_df = self.perturbedCheck()
                qual_assess_df = pd.merge(qual_assess_df, coverage_df, how='left', on='FASTQFILENAME')
        except AttributeError:
            self.logger.info('query_df or coverage_check_flag not present -- no coverage check')
        # format the qual_assess_df dataframe
        qual_assess_df = self.formatQualAssessDataFrame(qual_assess_df)
        return qual_assess_df

    @abc.abstractmethod
    def auditQualAssessDataFrame(self, qual_assess_df):
        """
            read in appropriate key/values from config file, add status and auto_audit columns to qual_asses_df
            :params qual_assess_df: the quality assessment dataframe with all appropriate columns for values in quality_assess thresholds in config file
            :returns: qual_assess_df with STATUS and AUTO_AUDIT columns added
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    @abc.abstractmethod
    def quantifyNonCodingRna(self, qual_assess_df):
        """
            quantify amount of non coding RNA in library. Note: this requires annotation of non coding RNA in annotation file
            :param qual_assess_df: qual_assess_df with at least FASTQFILENAME column (or some column used to identify each file uniquely)
            :returns qual_assess_df: noncoding_rna_df with columns 'FASTQFILENAME','TOTAL_rRNA', 'UNIQUE_rRNA', 'UNIQUE_tRNA_ncRNA' added.
            NOTE: FASTQFILENAME may be some other way of consistently identifying a sample uniquely
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    @abc.abstractmethod
    def formatQualAssessDataFrame(self, qual_assess_df):
        """
            A final step -- format columns in qual_assess_df for final writing
            :params qual_assess_df: a quality_assess_df with all columns
            :returns: qual_assess_df with appropriate formatting for writing out
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    def parseAlignmentLog(self, alignment_log_file_path):
        """
            parse the information on the alignment out of a novoalign log
            :param alignment_log_file_path: the filepath to a novoalign alignment log
            :returns: a dictionary of the parsed data of the input file
        """
        library_metadata_dict = {}
        alignment_regex_dict = {'LIBRARY_SIZE': r"(?<=Read Sequences:\s)\s*\d*",
                                'UNIQUE_ALIGNMENT': r"(?<=Unique Alignment:\s)\s*\d*",
                                'MULTI_MAP': r"(?<=Multi Mapped:\s)\s*\d*",
                                'NO_MAP': r"(?<=No Mapping Found:\s)\s*\d*",
                                'HOMOPOLY_FILTER': r"(?<=Homopolymer Filter:\s)\s*\d*",
                                'READ_LENGTH_FILTER': r"(?<=Read Length:\s)\s*\d*"}

        # open the log path
        alignment_file = open(alignment_log_file_path, 'r')
        alignment_file_text = alignment_file.read()
        # loop over alignment_regex dict and enter values extracted from alignment_file into alignment_metadata_dict
        for alignment_category, regex_pattern in alignment_regex_dict.items():
            # extract the value corresponding to the alignment_category regex (see alignment_regex_dict)
            try:
                extracted_value = int(re.findall(regex_pattern, alignment_file_text)[0])
            except ValueError:
                print('problem with file %s' % alignment_log_file_path)
            except IndexError:
                print('No %s in %s. Value set to 0' % (alignment_category, alignment_log_file_path))
                extracted_value = 0
            # check that the value is both an int and not 0
            if isinstance(extracted_value, int):
                library_metadata_dict.setdefault(alignment_category, extracted_value)
            else:
                print('cannot find %s in %s' % (alignment_category, alignment_log_file_path))

        # close the alignment_file and return
        alignment_file.close()
        return library_metadata_dict

    @abc.abstractmethod
    def parseGeneCount(self, htseq_counts_path):
        """
            count the gene counts that mapped either to genes (see COUNT_VARS at top of script for other features)
            :param htseq_counts_path: a path to a  _read_count.tsv file (htseq-counts output)
            :returns: a dictionary with the keys FEATURE_ALIGN_NOT_UNIQUE, TOO_LOW_AQUAL, AMBIGUOUS_FEATURE, NO_FEATURE, NOT_ALIGNED_TOTAL (at least, maybe more)
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    @abc.abstractmethod
    def uniqueAmbiguousProteinCodingCount(self, fastq_simplename):
        """
            intersect bed file with protein coding coords with an alignment file with the htseq annotations added (in this case, grep out __ambiguous)
            :params fastq_simplename: fastq filename minus any path and extention eg /path/to/my_reads_R1.fastq.gz would be my_reads_R1
            :returns: the number of reads (lines) aligning to protein coding coordinates
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    def totalrRNA(self, bam_path, rRNA_region, strandedness):
        """
            extract number of crypto aligning to rRNA
            :params bam_path: path to a bam file. requires that an index (with samtools index, novoindex, etc) be in same directory (will have .bai extension)
            :params rRNA_region: region of the rRNA in samtools view format -- chromosome:start_coord-stop_coord. eg CP022021.1:204523-235534
            :params strandedness: no or reverse, according to the library prep. This determines whether all
            :returns: total_rRNA is unique_rRNA + num_primary_alignment_rRNA (this is the number of primarily alignments from the multimap), AND unique_rRNA
        """
        # extract number of primary alignments of multimapped reads to rRNA. If reverse, only accept reads mapping in the positive direction
        # NOTE: THIS IS ONLY CORRECT IF THE rRNA IS ON THE FORWARD STRAND
        if strandedness == 'reverse':
            # -F 16 excludes reverse strand reads
            cmd_primary_multi_alignment_rRNA = 'samtools view -F 16 %s %s | grep ZS:Z:R | grep HI:i:1 | grep -v \"HI:i:1[[:digit:]]\" | wc -l' % (
                bam_path, rRNA_region)
        else:
            # grep -v excludes reads with a digit after the 1 (there is a prettier way to do this, im sure)
            cmd_primary_multi_alignment_rRNA = 'samtools view %s %s | grep ZS:Z:R | grep HI:i:1 | grep -v \"HI:i:1[[:digit:]]\" | wc -l' % (
                bam_path, rRNA_region)
        # as long as this is the first function called that needs an index, this will error check that samtools index has been run
        try:
            num_primary_alignment_rRNA = int(subprocess.getoutput(cmd_primary_multi_alignment_rRNA))
        except ValueError:
            sys.exit('You must first index the alignment files with samtools index')

        # extract number of unique alignments to rRNA
        # NOTE: THIS IS ONLY CORRECT IF THE rRNA IS ON THE FORWARD STRAND
        if strandedness == 'reverse':
            cmd_unique_rRNA = 'samtools view -F 16 %s %s | grep -v ZS:Z:R | wc -l' % (bam_path, rRNA_region)
        else:
            cmd_unique_rRNA = 'samtools view %s %s | grep -v ZS:Z:R | wc -l' % (bam_path, rRNA_region)
        unique_rRNA = int(subprocess.getoutput(cmd_unique_rRNA))

        # add for total rRNA
        total_rRNA = num_primary_alignment_rRNA + unique_rRNA

        return total_rRNA, unique_rRNA

    # TODO: CURRENTLY DOES NOT INCLUDE PRIMARY ALIGNMENT OF MULTI MAPS -- SHOULD IT?
    def totaltRNAncRNA(self, bam_path, trna_ncrna_annotation_gff, strandedness):
        """
            count reads that align to provided gff if the annotation overlaps at least 90% of the read uniquely
            :param bam_path: path to bam annotation file
            :param trna_ncrna_annotation_gff: a gff3 format with same chromosome/location notation as genome
            :param strandedness: strandedness of the library ('no' or 'reverse')
            :returns: number of reads fulfilling overlap criteria
        """
        # construct command based on strandedness of library
        if strandedness == 'reverse':
            # note: look up bedtools intersect --help for flags. grep -v returns all lines except those matching the pattern, in this case signifying multimaps
            bedtools_cmd = 'bedtools intersect -s -f .90 -a %s -b %s | samtools view | grep -v ZS:Z:R | wc -l' % (
                bam_path, trna_ncrna_annotation_gff)
        else:
            # no -s means this will count intersects regardless of strand
            bedtools_cmd = 'bedtools intersect -f .90 -a %s -b %s | samtools view | grep -v ZS:Z:R | wc -l' % (
                bam_path, trna_ncrna_annotation_gff)

        # extract unique_alignments to nc and t RNA
        unique_align_tRNA_ncRNA = int(subprocess.getoutput(bedtools_cmd))

        return unique_align_tRNA_ncRNA

    @abc.abstractmethod
    def calculateIntergenicCoverage(self, qual_assess_df):
        """
            calculate coverage of regions in genome between features
            Requires that the number of bases in the intergenic regions of the genome is present in OrganismData_config.ini
            in genome_files/<organism>
            TODO: BETTER NOTES HERE -- SEE CRYPTOQUALITYASSESSMENTOBJECT
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    @abc.abstractmethod
    def calculateExonicCoverage(self, bam_file, output_directory):
        """
            calculate coverage of exon regions. deposits file in output directory named utils.pathBaseName(bam_file)+'_exonic_coverage.csv'
            :param bam_file: path to bam file
            :param output_directory: path to output directory
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    def createCoverageBedFile(self, bam_file):
        """

        """
        raise NotImplementedError

    def calculatePercentFeatureCoverage(self, feature, genotype, annotation_path, bam_file, num_bases_in_region=None):
        """
            Calculate percent of given feature (regions summed, so all exon in gene, eg) of a gene (exon, CDS, etc) covered by 1 or more reads
            :param feature: annotation feature over which to take percentage, eg all exons in gene, or all CDS
            :param genotype: gene in annotation file
            :param annotation_path: path to annotation file
            :param bam_file: a sorted, indexed alignment file (.bam)
            :param num_bases_in_region: pass number of bases in the region directly, this will skip the step of calculating this number from the annotation file
            :returns: the fraction of bases in the (summed over the number of features in the gene) feature region covered by at least one read
        """
        if not num_bases_in_region:
            # extract number of bases in CDS of given gene. Credit: https://www.biostars.org/p/68283/#390427
            num_bases_in_region_cmd = "grep %s %s | grep %s | bedtools merge | awk -F\'\t\' \'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}\'" % (
                genotype, annotation_path, feature)
            num_bases_in_region = int(subprocess.getoutput(num_bases_in_region_cmd))
        # extract number of bases with depth != 0
        num_bases_depth_not_zero_cmd = "grep %s %s | grep %s | gff2bed | samtools depth -aa -Q 10 -b - %s | cut -f3 | grep -v 0 | wc -l" % (
            genotype, annotation_path, feature, bam_file)
        num_bases_in_cds_with_one_or_more_read = int(subprocess.getoutput(num_bases_depth_not_zero_cmd))

        return num_bases_in_cds_with_one_or_more_read / float(num_bases_in_region)

    def igvShot(self, bam_file_list, organism=None, igv_output_dir=None, gene=None, audited_qual_assess_df=None,
                manual_audit_flag=False):
        """
            first checks if alignment files are indexed; creates and submit sbatch script to index if not.
            creates a lookup_file.tsv of alignment_file(bam)_path \t bed_path \t igv_genome_path \t igv_output_dir
            sbatch file uses lookup file to input to this command:

            make_IGV_snapshots.py ${bam_file} -bin /opt/apps/igv/2.4.7/igv.jar -nf4 -r ${bed_file} -g ${igv_genome} -fig_format png -o ${igv_outout_dir}

            :param igv_output_dir: default None, and the output will be deposited in reports/igv_<datetime>. enter this to redirect output elsewhere
            :param bam_file_list: list of bam files to take igv shots of
            :param organism: default is none, but required if audited_qual_assess_df not passed
            :param gene: a specific gene, or list of genes, to take a browser shot of in all samples
            :param audited_qual_assess_df: a qual_assess_df with column AUTO_AUDIT and/or MANUAL_AUDIT filled with 0/1.
                                           If 1, a shot will be taken of the perturbed gene and wildtype in same cond. If no
                                           wt in same condition exists, the first wildtype in list will be used.
                                           May contain more than files in bam_file_list, but only those that can be found
                                           from the bam file will be used
            :param manual_audit_flag: default false. If true, column MANUAL_AUDIT will be used instead of AUTO_AUDIT
        """
        # ensure either gene or qudited_qual_assess_df is passed
        try:
            if not gene or audited_qual_assess_df:
                raise AttributeError
        except AttributeError:
            self.logger.critical(
                'igvShot cannot make lookup file without either a gene, gene list, or an audited_qual_assess_df')
            print('igvShot cannot make lookup file without either a gene, gene list, or an audited_qual_assess_df')
        # user must pass either gene/gene_list or a audited_qual_assess_df
        if gene is not None and audited_qual_assess_df is not None:
            raise ValueError('GeneAndQualAssessCannotBePassedTogether')
        # if gene is passed, the organism is required. In this case, the annotation file used will be the one in the annotation_file slot in OrganismData_config.ini in genome_files
        if gene and organism is None:
            raise AttributeError('IfGenePassedOrganismRequired')

        # set output_dir to default if not passed by user
        if not igv_output_dir:
            igv_output_dir = os.path.join(self.reports, 'igv_%s_%s' % (self.year_month_day, utils.hourMinuteSecond()))
            utils.mkdirp(igv_output_dir)

        # make sure that gene is a list (this allows user to pass multiple genes as a list, or single gene)
        if not isinstance(gene, list):
            gene = [gene]

        # check to see if index files are present. If not, create a list to pass to indexBamFile
        bam_files_to_index = []
        for bam_file in bam_file_list:
            if not os.path.isfile(bam_file + '.bai'):
                bam_files_to_index.append(bam_file)

        # if there are unindexed bam files, notify user and create batchscript to index them
        try:
            if not len(bam_files_to_index) == 0:
                raise AttributeError('NotAllBamFilesIndexed')
        except AttributeError:
            # create sbatch lookup and script
            self.indexBamFileBatchScript(bam_files_to_index)
            print('Some bam files not indexed. Creating batchscript to index these files.\n'
                  'Submit the batchscript to index the bamfiles and then re-run igvShot')

        # create lookup file for sbatch
        if audited_qual_assess_df:
            lookup_file_path = self.createIgvLookupFromDataframe(bam_file_list, audited_qual_assess_df,
                                                                 manual_audit_flag)
        else:
            lookup_file_path = self.createIgvLookupFromGeneList(organism, bam_file_list,
                                                                gene)  # note 'gene' may be a list, and even if it is not passed as one, will be turned into a list of length 1

        # create sbatch script for igv shots
        self.createIgvBatchscript(lookup_file_path, igv_output_dir)

    def extractLog2cpm(self, gene, fastq_simple_name, log2cpm_csv_path):
        """
           extract log2cpm
           :param gene: log2cpm of gene you wish to extract
           :param fastq_simple_name: name of the fastq_file, no ext, no path
           :param log2cpm_csv_path: created by raw_counts.py --> log2_cpm.R. index in gene_name, columns are fastq_simple_name_read_count.tsv
           :returns: log2cpm of a gene in a given library
        """
        # read in log2cpm dataframe
        try:
            log2cpm_df = utils.readInDataframe(log2cpm_csv_path)
            log2cpm_df = log2cpm_df.set_index('gene_id')
        except FileNotFoundError:
            error_msg = 'path to log2_cpm not valid: %s' %log2cpm_csv_path
            self.logger.critical(error_msg)
            print(error_msg)

        # create column name corresponding to log2cpm_df from fastq_simple_name
        column_name = fastq_simple_name + '_read_count.tsv'
        # check that it is actually in the log2cpm_df columns
        try:
            if column_name not in log2cpm_df.columns:
                raise AttributeError('ColumnNameNotInLog2CpmSheet')
        except AttributeError:
            error_msg = '%s not in log2cpm sheet %s' %(column_name, log2cpm_csv_path)
            self.logger.critical(error_msg)
            print(error_msg)
        else:
            # check that the gene is in the gene_id column
            try:
                if gene not in log2cpm_df.index:
                    gene = gene.replace('CNAG', 'CKF44')
                    if gene not in log2cpm_df.index:
                        raise AttributeError('GeneNotInLog2cpmSheet')
            except AttributeError:
                error_msg = '%s not in log2cpm sheet %s' %(gene, log2cpm_csv_path)
                self.logger.critical(error_msg)
                print(error_msg)

            else:
                # extract log2cpm and return
                return log2cpm_df.loc[gene, column_name]



    def indexBamFileBathScript(self, bam_files_to_index):
        """

        """
        raise NotImplementedError('AbstractMethodMustBeOverwritten')

    def createIgvLookupFromDataframe(self, bam_list, audited_qual_assess_df, manual_audit_flag):
        """

        :returns: path to the lookup file for createIgvBatchscript
        """
        # check inputs
        # set default wildtype as first wildtype in list (use this for wildtype if no wildtype in same cond/timepoint present

        # decompose bit status if manual_audit_flag is false, only take screen shots of perturbation, overexpression fail + drug markers and wildtype,
        # wildtypes with drug markers if drug marker wt fail, and perturbation + drug markers for perturbation drug marker fail
        # if wildtype of same condition, timepoint, etc doesn't work, use default wildtype
        # use createIgvLookupFromGeneList and conconcatenate?

        raise NotImplementedError('AbstractMethodMustBeOverwritten')

    def createIgvLookupFromGeneList(self, organism, bam_list, gene_list, gene_offset=500):
        """

        :returns: path to the lookup file for createIgvBatchscript
        """
        organism_genome_files = os.path.join(self.genome_files, organism)
        if not os.path.isdir(organism_genome_files):
            self.logger.critical('%s not a in genome_files. Make sure genome_files exists in rnaseq_pipeline, '
                                 'and that the organism is passed as a correctly formatted subdirectory of genome_files' % organism_genome_files)
            raise NotADirectoryError('OrganismDirectoryNotFoundInGenomeFiles')
        organism_config_file = os.path.join(organism_genome_files, 'OrganismData_config.ini')
        if not os.path.isfile(organism_config_file):
            self.logger.critical(
                '%s not found -- check genome_files/organism subdir. Possibly delete genome_files and let the script re-download (make sure it is accessible)' % organism_config_file)
            raise FileNotFoundError('ConfigFileNotFound')
        # read in organism config file
        organism_config = configparser.ConfigParser()
        organism_config.read(organism_config_file)
        organism_config_dict = organism_config['OrganismData']
        # set annotation file
        annotation_file = os.path.join(self.genome_files, organism, organism_config_dict['annotation_file'])
        if not os.path.isfile(annotation_file):
            self.logger.critical('the annotation file %s not valid filepath' % annotation_file)
            raise FileNotFoundError('AnnotationFileNotFound')
        # set igv genome
        try:
            igv_genome = os.path.join(self.genome_files, organism, organism_config_dict['igv_genome'])
        except KeyError:  # this is for KN99
            igv_genome = os.path.join(self.genome_files, organism, organism_config_dict['igv_stranded_genome'])
        if not os.path.isfile(igv_genome):
            self.logger.critical('igv_genome path %s not valid' % igv_genome)
            raise FileNotFoundError('IgvGenomeNotFound')
        # make directory to store bedfiles
        bed_file_dir_path = os.path.join(self.job_scripts, 'igv_%s' % self.year_month_day)
        utils.mkdirp(bed_file_dir_path)
        # make lookup_file_path
        igv_lookup_file_path = os.path.join(bed_file_dir_path, 'igv_lookup_file.txt')

        bed_entry_dict = {}
        for gene in gene_list:
            # get first column corresponding to given gene, take uniq value as chromosome
            extract_chr_cmd = 'grep %s %s | cut -f1 | uniq' % (gene, annotation_file)
            chromosome_identifier = subprocess.getoutput(extract_chr_cmd)

            # get the list of all start coordinates associated with a feature, sort, and take smallest
            extract_start_coord_cmd = 'grep %s %s | cut -f4 | sort -n | head -1' % (gene, annotation_file)
            start_coord = int(subprocess.getoutput(extract_start_coord_cmd))
            start_coord_with_offest = start_coord - gene_offset  # add offset to widen window on igv_viewer
            if start_coord_with_offest < 1:
                start_coord_with_offest = 1
            # sort stop coordinates of all features of given gene, take largest
            extract_stop_coord_cmd = 'grep %s %s | cut -f5 | sort -n | tail -1' % (gene, annotation_file)
            stop_coord = int(subprocess.getoutput(extract_stop_coord_cmd))
            stop_coord_with_offset = stop_coord + gene_offset  # add offset to widen window on igv_viewer

            # enter {gene: bed_line} to dict -- not the final \t allows to add the bam_file_simplename later
            bed_line = '%s\t%s\t%s\t' % (chromosome_identifier, start_coord_with_offest, stop_coord_with_offset)
            bed_entry_dict.setdefault(gene, bed_line)

        for bam_file in bam_list:
            bam_simple_name = utils.pathBaseName(bam_file)
            for gene in bed_entry_dict:
                bed_file_path = os.path.join(bed_file_dir_path, bam_simple_name + '_%s.bed' % gene)
                with open(bed_file_path, 'w') as bed_file:
                    bed_file.write(bed_entry_dict[gene] + '%s' % bam_simple_name)
                with open(igv_lookup_file_path, 'a') as igv_lookup_file:
                    new_lookup_line = '%s\t%s\t%s\n' % (bam_file, bed_file_path, igv_genome)
                    igv_lookup_file.write(new_lookup_line)

        return igv_lookup_file_path

    def createIgvBatchscript(self, lookup_file_path, output_dir):
        """

        """
        # error check input
        if not os.path.isfile(lookup_file_path):
            self.logger.critical('lookup file %s does not exist' % lookup_file_path)
            raise FileNotFoundError('LookupFileDoesNotExist')
        if not os.path.isdir(output_dir):
            self.logger.critical('output dir %s does not exist' % output_dir)
            raise NotADirectoryError('OutputDirectoryDoesNotExist')

        # write sbatch job. see https://htcfdocs.readthedocs.io/en/latest/runningjobs/
        line_count_cmd = 'cat %s | wc -l' % lookup_file_path
        line_count = int(subprocess.getoutput(line_count_cmd))
        sbatch_array_line = "--array=1-{}%{}".format(line_count, min(line_count, 20))
        job = '#!/bin/bash\n\n' \
              '#SBATCH -N 1\n' \
              '#SBATCH --mem=20G\n' \
              '#SBATCH %s\n' \
              '#SBATCH -o %s/igv_snapshot_%s_%s.out\n' \
              '#SBATCH -J igv_snapshot\n\n' \
              'ml rnaseq_pipeline java\n\n' \
              'read bam_file bed_file igv_genome < <(sed -n ${SLURM_ARRAY_TASK_ID}p %s )\n\n' \
              'make_IGV_snapshots.py ${bam_file} -bin /opt/apps/igv/2.4.7/igv.jar -nf4 -r ${bed_file} -g ${igv_genome} -fig_format png -o %s\n' \
              % (sbatch_array_line, self.sbatch_log, self.year_month_day, utils.hourMinuteSecond(), lookup_file_path,
                 output_dir)

        # write to file
        # check that a job_scripts directory for igv for today has been created, make one if no
        job_script_dir_path = os.path.join(self.job_scripts, 'igv_%s' % self.year_month_day)
        if not os.path.isdir(os.path.join(job_script_dir_path)):
            utils.mkdirp(job_script_dir_path)

        igv_job_script_path = os.path.join(job_script_dir_path, 'igv_%s_%s.sbatch' % (self.year_month_day, utils.hourMinuteSecond()))
        with open(igv_job_script_path, 'w') as igv_job_script_file:
            igv_job_script_file.write(job)

        # error check igv_job_script_path
        if not os.path.isfile(igv_job_script_path):
            self.loger.critical('Failed to write job %s to file %s' % (job, igv_job_script_path))
            raise FileNotFoundError('IgvJobScriptNotSuccessfullyWritten')
        else:
            print('Submitting igv batchscript %s' % igv_job_script_path)
            batchscript_submit_cmd = 'sbatch %s' % igv_job_script_path
            utils.executeSubProcess(batchscript_submit_cmd)

    def qortsPlots(self):
        """
            create the qorts multiplots
        """
        try:  # TODO: clean this up
            if not os.path.isfile(self.query_path):
                raise FileExistsError('QueryPathNotValid')
            else:
                if not hasattr(self, 'standardized_database_df'):
                    query_df = utils.readInDataframe(self.query_path)
                    self.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(query_df)
        except FileExistsError or AttributeError:
            print('Query Path not valid')

        try:
            if not os.path.isdir(self.qorts_output):
                raise NotADirectoryError('QortsOutputDoesNotExist')
        except NotADirectoryError or AttributeError:
            print('qorts_output is either not passed as an argument, or does not exist.')

        try:
            if not hasattr(self, 'organism'):
                raise AttributeError('NoOrganismPassed')
        except AttributeError:
            print('The QualityAssessmentObject must have an attribute organism')

        try:
            if not os.path.isfile(self.annotation_file):
                raise FileExistsError('AnnotationFilePathNotValid')
        except (FileExistsError, AttributeError):
            organism_object = OrganismData(organism=self.organism)

        try:
            if not hasattr(self, 'bam_file_list'):
                raise AttributeError('NoBamFileList')
        except AttributeError:
            print('You must supply a list of bamfiles')


        # sort list of countfilenames pre and post 2015
        self.standardized_database_df['LIBRARYDATE'] = pd.to_datetime(self.standardized_database_df['LIBRARYDATE'])
        pre_2015_df = self.standardized_database_df[self.standardized_database_df['LIBRARYDATE'] <= '2015-01-01']
        pre_2015_sample_list = list(pre_2015_df['FASTQFILENAME'])
        post_2015_df = self.standardized_database_df[self.standardized_database_df['LIBRARYDATE'] > '2015-01-01']
        post_2015_sample_list = list(post_2015_df['FASTQFILENAME'])

        cmd_list = []

        for fastqfile_simplename in post_2015_sample_list:
            bamfilename_path = [bam_file for bam_file in self.bam_file_list if fastqfile_simplename in bam_file][0]
            try:
                if bamfilename_path not in self.bam_file_list:
                    raise FileNotFoundError('QuerySampleNotInAlignCountDirectory')
            except FileNotFoundError:
                return ('A sample in the query could not be located in the directory with the alignment files.\n'
                        'Make sure the query provided corresponds to the align_count_path directory')

            output_subdir = os.path.join(self.qorts_output, fastqfile_simplename + '_qorts')
            utils.mkdirp(output_subdir)

            try:
                if not os.path.exists(bamfilename_path):
                    raise FileExistsError('BamfilePathNotValid')
            except FileExistsError:
                print('path from align_count_path to bamfile does not exist')
            if not hasattr(self, 'annotation_file'):
                annotation_file = organism_object.annotation_file
            else:
                annotation_file = self.annotation_file

            try:
                if not os.path.exists(bamfilename_path):
                    raise FileExistsError('BamfilePathNotValid')
            except FileExistsError:
                print('path from align_count_path to bamfile does not exist')
            qorts_cmd = 'java -Xmx1G -jar /opt/apps/labs/mblab/software/hartleys-QoRTs-099881f/scripts/QoRTs.jar QC --singleEnded --stranded --keepMultiMapped --generatePlots %s %s %s\n' % (
                bamfilename_path, annotation_file, output_subdir)
            cmd_list.append(qorts_cmd)

        for fastqfile_simplename in pre_2015_sample_list:
            bamfilename_path = [bam_file for bam_file in self.bam_file_list if fastqfile_simplename in bam_file][0]
            try:
                if bamfilename_path not in self.bam_file_list:
                    raise FileNotFoundError('QuerySampleNotInAlignCountDirectory')
            except FileNotFoundError:
                print('A sample in the query could not be located in the directory with the alignment files.\n'
                      'Make sure the query provided corresponds to the align_count_path directory')

            output_subdir = os.path.join(self.qorts_output, fastqfile_simplename + '_qorts')
            utils.mkdirp(output_subdir)

            try:
                if not os.path.exists(bamfilename_path):
                    raise FileExistsError('BamfilePathNotValid')
            except FileExistsError:
                print('path from align_count_path to bamfile does not exist')
            if not hasattr(self, 'annotation_file'):
                if hasattr(organism_object, 'annotation_file_no_strand'):
                    annotation_file = organism_object.annotation_file_no_strand
                else:
                    annotation_file = organism_object.annotation_file
            else:
                annotation_file = self.annotation_file

            qorts_cmd = 'java -Xmx1G -jar /opt/apps/labs/mblab/software/hartleys-QoRTs-099881f/scripts/QoRTs.jar QC --singleEnded --keepMultiMapped --generatePlots %s %s %s\n' % (
                bamfilename_path, annotation_file, output_subdir)
            cmd_list.append(qorts_cmd)

        with open(self.sbatch_script, 'w') as sbatch_file:
            sbatch_file.write('\n')
            sbatch_file.write('\n')

            for cmd in cmd_list:
                sbatch_file.write(cmd)
                sbatch_file.write('\n')
