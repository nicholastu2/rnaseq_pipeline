import os
import subprocess
import re
import pandas as pd
from glob import glob
import sys
from rnaseq_tools import utils
from rnaseq_tools.StandardDataObject import StandardData
from rnaseq_tools.DatabaseObject import DatabaseObject
import abc

# turn off SettingWithCopyWarning in pandas
pd.options.mode.chained_assignment = None

class QualityAssessmentObject(StandardData):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = ['bam_file_list', 'count_file_list', 'novoalign_log_list', 'coverage_check_flag',
                                         'query_path', 'standardized_database_df', 'qual_assess_dir_path']
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
            library_metadata_dict['AMBIGUOUS_UNIQUE_PROTEIN_CODING_READS'] = self.uniqueAmbiguousProteinCodingCount(fastq_basename)
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
        raise NotImplementedError

    @abc.abstractmethod
    def quantifyNonCodingRna(self, qual_assess_df):
        """
            quantify amount of non coding RNA in library. Note: this requires annotation of non coding RNA in annotation file
            :param qual_assess_df: qual_assess_df with at least FASTQFILENAME column (or some column used to identify each file uniquely)
            :returns qual_assess_df: noncoding_rna_df with columns 'FASTQFILENAME','TOTAL_rRNA', 'UNIQUE_rRNA', 'UNIQUE_tRNA_ncRNA' added.
            NOTE: FASTQFILENAME may be some other way of consistently identifying a sample uniquely
        """
        raise NotImplementedError

    @abc.abstractmethod
    def formatQualAssessDataFrame(self, qual_assess_df):
        """
            A final step -- format columns in qual_assess_df for final writing
            :params qual_assess_df: a quality_assess_df with all columns
            :returns: qual_assess_df with appropriate formatting for writing out
        """
        raise NotImplementedError

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
        raise NotImplementedError

    @abc.abstractmethod
    def uniqueAmbiguousProteinCodingCount(self, fastq_simplename):
        """
            intersect bed file with protein coding coords with an alignment file with the htseq annotations added (in this case, grep out __ambiguous)
            :params fastq_simplename: fastq filename minus any path and extention eg /path/to/my_reads_R1.fastq.gz would be my_reads_R1
            :returns: the number of reads (lines) aligning to protein coding coordinates
        """
        raise NotImplementedError

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
            cmd_primary_multi_alignment_rRNA = 'samtools view -F 16 %s %s | grep ZS:Z:R | grep HI:i:1 | grep -v \"HI:i:1[[:digit:]]\" | wc -l' % (bam_path, rRNA_region)
        else:
            # grep -v excludes reads with a digit after the 1 (there is a prettier way to do this, im sure)
            cmd_primary_multi_alignment_rRNA = 'samtools view %s %s | grep ZS:Z:R | grep HI:i:1 | grep -v \"HI:i:1[[:digit:]]\" | wc -l' % (bam_path, rRNA_region)
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
            bedtools_cmd = 'bedtools intersect -s -f .90 -a %s -b %s | samtools view | grep -v ZS:Z:R | wc -l' % (bam_path, trna_ncrna_annotation_gff)
        else:
            # no -s means this will count intersects regardless of strand
            bedtools_cmd = 'bedtools intersect -f .90 -a %s -b %s | samtools view | grep -v ZS:Z:R | wc -l' % (bam_path, trna_ncrna_annotation_gff)

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
        raise NotImplementedError

    @abc.abstractmethod
    def calculateExonicCoverage(self, bam_file, output_directory):
        """
            calculate coverage of exon regions. deposits file in output directory named utils.pathBaseName(bam_file)+'_exonic_coverage.csv'
            :param bam_file: path to bam file
            :param output_directory: path to output directory
        """

    def calculatePercentFeatureCoverage(self, feature, genotype, annotation_path, bam_file, num_bases_in_region = None):
        """
            Calculate percent of given feature (regions summed, so all exon in gene, eg) of a gene (exon, CDS, etc) covered by 1 or more reads
            :param feature: annotation feature over which to take percentage, eg all exons in gene, or all CDS
            :param genotype: gene in annotation file
            :param bam_file: a sorted, indexed alignment file (.bam)
            :param num_bases_in_region: pass number of bases in the region directly, this will skip the step of calculating this number from the annotation file
            :returns: the fraction of bases in the (summed over the number of features in the gene) feature region covered by at least one read
        """
        if not num_bases_in_region:
            # extract number of bases in CDS of given gene. Credit: https://www.biostars.org/p/68283/#390427
            num_bases_in_region_cmd = "grep %s %s | grep %s | bedtools merge | awk -F\'\t\' \'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}\'" % (genotype, annotation_path, feature)
            num_bases_in_region = int(subprocess.getoutput(num_bases_in_region_cmd))
        # extract number of bases with depth != 0
        num_bases_depth_not_zero_cmd = "grep %s %s | grep %s | gff2bed | samtools depth -a -b - %s | cut -f3 | grep -v 0 | wc -l" %(genotype, annotation_path, feature, bam_file)
        num_bases_in_cds_with_one_or_more_read = int(subprocess.getoutput(num_bases_depth_not_zero_cmd))

        return num_bases_in_cds_with_one_or_more_read / float(num_bases_in_region)

    def igvShot(self):
        """
            create IgvObject and create browser shots
        """
        raise NotImplementedError

    def qortsPlots(self):
        """
            create the qorts multiplots
        """
        try:
            if not os.path.isfile(self.align_count_path):
                raise NotADirectoryError('AlignCountPathNotValid')
        except NotADirectoryError or AttributeError:
            print(
                'Align count path is either not set or not valid. align_count_path should be pointed toward a directory\n'
                'with alignment files, typically the output of align_count.py or create_experiment.py')

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
            if not os.path.isfile(self.annotation_file):
                raise FileExistsError('AnnotationFilePathNotValid')
        except FileExistsError or AttributeError:
            print('You must pass a correct path to an annotation file (see genome_files)')

        sorted_bamfile_list = glob(os.path.join(self.align_count_path, '*_sorted_aligned_reads.bam'))
        sorted_bamfile_list = [os.path.basename(x) for x in sorted_bamfile_list]

        # sort list of countfilenames pre and post 2015
        self.standardized_database_df['LIBRARYDATE'] = pd.to_datetime(self.standardized_database_df['LIBRARYDATE'])
        pre_2015_df = self.standardized_database_df[self.standardized_database_df['LIBRARYDATE'] <= '2015-01-01']
        pre_2015_countfilename_list = list(pre_2015_df['COUNTFILENAME'])
        post_2015_df = self.standardized_database_df[self.standardized_database_df['LIBRARYDATE'] > '2015-01-01']
        post_2015_countfilename_list = list(post_2015_df['COUNTFILENAME'])

        cmd_list = []

        for countfilename in post_2015_countfilename_list:
            bamfilename = countfilename.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
            try:
                if not bamfilename in sorted_bamfile_list:
                    raise FileNotFoundError('QuerySampleNotInAlignCountDirectory')
            except FileNotFoundError:
                return ('A sample in the query could not be located in the directory with the alignment files.\n'
                        'Make sure the query provided corresponds to the align_count_path directory')

            output_subdir = os.path.join(self.qorts_output, '[' + bamfilename + ']' + '_qorts')
            utils.mkdirp(output_subdir)
            bamfilename_path = os.path.join(self.align_count_path, bamfilename)
            try:
                if not os.path.exists(bamfilename_path):
                    raise FileExistsError('BamfilePathNotValid')
            except FileExistsError:
                print('path from align_count_path to bamfile does not exist')
            qorts_cmd = 'java -Xmx1G -jar /opt/apps/labs/mblab/software/hartleys-QoRTs-099881f/scripts/QoRTs.jar QC --singleEnded --stranded --keepMultiMapped --generatePlots %s %s %s\n' % (
                bamfilename_path, self.annotation_file, output_subdir)
            cmd_list.append(qorts_cmd)

        for countfilename in pre_2015_countfilename_list:
            bamfilename = countfilename.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
            try:
                if not bamfilename in sorted_bamfile_list:
                    raise FileNotFoundError('QuerySampleNotInAlignCountDirectory')
            except FileNotFoundError:
                print('A sample in the query could not be located in the directory with the alignment files.\n'
                      'Make sure the query provided corresponds to the align_count_path directory')

            output_subdir = os.path.join(self.qorts_output, '[' + bamfilename + ']' + '_qorts')
            utils.mkdirp(output_subdir)
            bamfilename_path = os.path.join(self.align_count_path, bamfilename)
            try:
                if not os.path.exists(bamfilename_path):
                    raise FileExistsError('BamfilePathNotValid')
            except FileExistsError:
                print('path from align_count_path to bamfile does not exist')
            qorts_cmd = 'java -Xmx1G -jar /opt/apps/labs/mblab/software/hartleys-QoRTs-099881f/scripts/QoRTs.jar QC --singleEnded --keepMultiMapped --generatePlots %s %s %s\n' % (
                bamfilename_path, self.annotation_file, output_subdir)
            cmd_list.append(qorts_cmd)

        with open(self.sbatch_script, 'w') as sbatch_file:
            sbatch_file.write('\n')
            sbatch_file.write('\n')

            for cmd in cmd_list:
                sbatch_file.write(cmd)
                sbatch_file.write('\n')