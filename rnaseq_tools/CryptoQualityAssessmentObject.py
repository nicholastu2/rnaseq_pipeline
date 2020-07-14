import os
import subprocess
import re
import pandas as pd
import configparser
from glob import glob
import sys
from rnaseq_tools import utils
from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject
from rnaseq_tools.DatabaseObject import DatabaseObject
import numpy as np

# turn off SettingWithCopyWarning in pandas
pd.options.mode.chained_assignment = None


class CryptoQualityAssessmentObject(QualityAssessmentObject):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = ['ko_gene_list', 'overexpress_gene_list', ]
        # This is a method of adding expected attributes to StandardData from StandardData children
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(CryptoQualityAssessmentObject, self).__init__(self._add_expected_attributes, **kwargs)
        # set standardDirectory structure
        self.standardDirectoryStructure()
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'QualityAssessmentObject'
        try:
            self.ko_gene_list = kwargs['ko_gene_list']
        except KeyError:
            self.ko_gene_list = []  # eg ['CNAG_01020' , ['CNAG_39392','CNAG_48382'], 'CNAG_23421'] where the center item is a double ko
        try:
            self.overexpress_gene_list = kwargs['overexpress_gene_list']
        except KeyError:
            self.overexpress_gene_list = []  # expecting no nested lists in this

    def auditQualAssess1(self, qual_assess_df):
        """
            use rnaseq_pipeline/config/quality_assess_config.ini entries to add status, auto_audit columns
            :params qual_assess_df: a quality_assess_df created by qual_assess_1 or a path to one
            :returns: qual_assess_df with added status and auto_audit columns
        """
        if os.path.isfile(qual_assess_df):
            qual_assess_df = utils.readInDataframe(qual_assess_df)

        # extract threshold/status from config file
        qual_assess_config = configparser.ConfigParser()
        qual_assess_config.read(self.config_file)
        qual_assess_1_dict = qual_assess_config['CryptoQualityAssessOne']

        # extract thresholds and bit status
        protein_coding_total_threshold = int(qual_assess_1_dict['PROTEIN_CODING_TOTAL_THRESHOLD'])
        protein_coding_total_bit_status = int(qual_assess_1_dict['PROTEIN_CODING_TOTAL_STATUS'])
        coverage_threshold = float(qual_assess_1_dict['COVERAGE_THRESHOLD'])
        coverage_bit_status = int(qual_assess_1_dict['COVERAGE_STATUS'])

        status_column_list = []
        for index, row in qual_assess_df.iterrows():
            # extract quality_assessment_metrics
            library_protein_coding_total = int(row['PROTEIN_CODING_TOTAL'])
            if row['GENOTYPE_1_COVERAGE'] is not None:
                library_genotype_1_coverage = float(row['GENOTYPE_1_COVERAGE'])
            else:
                library_genotype_1_coverage = -1
            if row['GENOTYPE_2_COVERAGE'] is not None:
                library_genotype_2_coverage = float(row['GENOTYPE_2_COVERAGE'])
            else:
                library_genotype_2_coverage = -1

            # set status_total to 0
            status_total = 0
            if library_protein_coding_total < protein_coding_total_threshold:
                status_total += protein_coding_total_bit_status
            if library_genotype_1_coverage > coverage_threshold or library_genotype_2_coverage > coverage_threshold:
                status_total += coverage_bit_status
            status_column_list.append(status_total)
            # TODO: ADD OVEREXPRESSION AND MARKER

        qual_assess_df['STATUS'] = status_column_list
        qual_assess_df['AUTO_AUDIT'] = np.where(qual_assess_df.STATUS > 0, 1, 0)

        return qual_assess_df

    def quantifyNonCodingRna(self, qual_assess_df):
        """

        """
        num_reads_to_ncRNA_dict = {}
        # set threshold to determine strandedness. note that there is an email from holly to yiming mentioning 10.25.2015. That is the best record we have, if htis message remains
        strandedness_date_threshold = pd.to_datetime('10.01.2015')
        kn99_tRNA_ncRNA_annotations = os.path.join(self.genome_files, 'KN99', 'ncRNA_tRNA_no_rRNA.gff')
        if hasattr(self, 'query_df'):
            for index, row in qual_assess_df.iterrows():
                genotype = list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME'] + '.fastq.gz')]['genotype'])[0]
                if genotype.startswith('CNAG'):
                    # extract fastq_filename without any preceeding path or file extension
                    fastq_simple_name = utils.pathBaseName(row['FASTQFILENAME'])
                    print('...evaluating ncRNA in %s' %fastq_simple_name)
                    # use this to extract bam_path
                    try:
                        bam_path = [bam_file for bam_file in self.bam_file_list if fastq_simple_name in bam_file][0]
                    except IndexError:
                        self.logger.info('%s not in bam_file_list' % fastq_simple_name)
                        continue
                    try:
                        if not os.path.isfile(bam_path):
                            raise FileNotFoundError
                    except FileNotFoundError:
                        self.logger.error('bam file not found %s' % bam_path)
                        print('bam file not found: %s' % bam_path)
                    libraryDate = list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME'] + '.fastq.gz')]['libraryDate'])[0]
                    row_date_time = pd.to_datetime(libraryDate)
                    strandedness = 'no' if row_date_time < strandedness_date_threshold else 'reverse'
                    total_rRNA, unique_rRNA = self.totalrRNA(bam_path, 'CP022322.1:272773-283180', strandedness)
                    unique_tRNA_ncRNA = self.totaltRNAncRNA(bam_path, kn99_tRNA_ncRNA_annotations, strandedness)
                    num_reads_to_ncRNA_dict.setdefault(fastq_simple_name, {'total_rRNA': total_rRNA, 'unique_rRNA': unique_rRNA, 'total_tRNA_ncRNA': unique_tRNA_ncRNA})

        # create dataframe from num_reads_to_ncRNA_dict
        ncRNA_df = pd.DataFrame.from_dict(num_reads_to_ncRNA_dict, orient='index').reset_index()
        # format column headers
        ncRNA_df.columns = ['FASTQFILENAME', 'TOTAL_rRNA', 'UNIQUE_rRNA', 'UNIQUE_tRNA_ncRNA']

        return ncRNA_df

    def formatQualAssess1DataFrame(self, qual_assess_df):
        """

        """
        # EFFECTIVE_LIBRARY_SIZE is LIBRARY_SIZE - (total_rRNA + unique_tRNA_ncRNA)
        qual_assess_df['EFFECTIVE_LIBRARY_SIZE'] = qual_assess_df['LIBRARY_SIZE'].astype('float') - (qual_assess_df['TOTAL_rRNA'] + qual_assess_df['UNIQUE_tRNA_ncRNA'])

        # EFFECTIVE_UNIQUE_ALIGNMENT is number of unique reads minus those unique reads mapping to rRNA and nc + t RNA
        qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT'] = qual_assess_df['UNIQUE_ALIGNMENT'] - (qual_assess_df['UNIQUE_rRNA'] +
                                                                                             qual_assess_df['UNIQUE_tRNA_ncRNA'] +
                                                                                             qual_assess_df['TOO_LOW_AQUAL'])

        # present the following categories as fraction of library_size
        # percent of library made up of rRNA (recall total_rRNA is unique + primary alignment since reads multimap in two spots, both seemingly rRNA)
        qual_assess_df['rRNA_PERCENT'] = qual_assess_df['TOTAL_rRNA'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
        qual_assess_df['nctRNA_PERCENT'] = qual_assess_df['UNIQUE_tRNA_ncRNA'] / qual_assess_df['LIBRARY_SIZE'].astype(float)
        qual_assess_df['nctrRNA_PERCENT'] = (qual_assess_df['TOTAL_rRNA'] + qual_assess_df['UNIQUE_tRNA_ncRNA']) / qual_assess_df['LIBRARY_SIZE'].astype(float)
        qual_assess_df['MULTI_MAP_PERCENT'] = qual_assess_df['MULTI_MAP'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
        qual_assess_df['NO_MAP_PERCENT'] = qual_assess_df['NO_MAP'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
        qual_assess_df['HOMOPOLY_FILTER_PERCENT'] = qual_assess_df['HOMOPOLY_FILTER'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
        qual_assess_df['READ_LENGTH_FILTER_PERCENT'] = qual_assess_df['READ_LENGTH_FILTER'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
        # htseq output not_aligned_total_percent is no_map + homopoly_filter + read_length filter. present as fraction of library_size
        qual_assess_df['NOT_ALIGNED_TOTAL_PERCENT'] = qual_assess_df['NOT_ALIGNED_TOTAL'] / qual_assess_df['LIBRARY_SIZE'].astype('float')

        # PROTEIN_CODING_TOTAL (formerly with_feature) is the number of reads mapping to a protein coding gene by htseq plus the unique ambiguous reads mapping to exon portiosn of overlapping protein coding reads
        qual_assess_df['PROTEIN_CODING_TOTAL'] = qual_assess_df['PROTEIN_CODING_COUNTED'] + qual_assess_df['AMBIGUOUS_UNIQUE_PROTEIN_CODING_READS']

        # protein_coding_total as percent of effective unique alignment
        qual_assess_df['PROTEIN_CODING_TOTAL_PERCENT'] = qual_assess_df['PROTEIN_CODING_TOTAL'] / qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT'].astype('float')
        # protein_coding_counted as percent of unique alignment
        qual_assess_df['PROTEIN_CODING_COUNTED_PERCENT'] = qual_assess_df['PROTEIN_CODING_COUNTED'] / qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT'].astype('float')

        # present the following as fraction of (total) unique_alignment  ### TODO: SHOULD THIS BE OF EFFECTIVE_UNIQUE_ALIGNMENT? PROBABLY YES: NOTE THE GRAPHS IN CURRENT OLD CRYPTO SHEETS ARE OF UNIQUE_ALIGNMENT
        qual_assess_df['NO_FEATURE_PERCENT'] = qual_assess_df['NO_FEATURE'] / qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT'].astype(float)
        qual_assess_df['AMBIGUOUS_FEATURE_PERCENT'] = qual_assess_df['AMBIGUOUS_FEATURE'] / qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT'].astype(float)
        qual_assess_df['TOO_LOW_AQUAL_PERCENT'] = qual_assess_df['TOO_LOW_AQUAL'] / qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT'].astype(float)

        # present EFFECTIVE_UNIQUE_ALIGNMENT as percent of library size (make sure this is the last step
        qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT_PERCENT'] = qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT'] / qual_assess_df['LIBRARY_SIZE'].astype(float)

        return qual_assess_df

    def parseGeneCount(self, htseq_counts_path):
        """
            NOTE: SPECIFICALLY SET UP FOR CRYPTO
            count the gene counts that mapped either to genes (see COUNT_VARS at top of script for other features)
            :param htseq_counts_path: a path to a  _read_count.tsv file (htseq-counts output)
            :returns: a dictionary with the keys FEATURE_ALIGN_NOT_UNIQUE, TOO_LOW_AQUAL, AMBIGUOUS_FEATURE, NO_FEATURE, NOT_ALIGNED_TOTAL
        """

        library_metadata_dict = {}
        # TODO: error checking on keys
        htseq_file = open(htseq_counts_path, 'r')
        htseq_file_reversed = reversed(htseq_file.readlines())

        crypto_protein_coding_count = 0
        line = next(htseq_file_reversed)
        try:
            while True:
                line_strip_split = line.strip().split('\t')
                if line.startswith('CKF44'):
                    # split the line, take the entry in the second column, which is the gene count, and add to crypto_protein_coding_effective_count
                    gene_count = int(line_strip_split[1])
                    crypto_protein_coding_count += gene_count
                if not (line.startswith('CNAG') or line.startswith('CKF44')):
                    # strip newchar, split on tab
                    line = line.strip().split('\t')
                    # extract the category of metadata count (eg __alignment_not_unique --> ALIGNMENT_NOT_UNIQUE)
                    htseq_count_metadata_category = line_strip_split[0][2:].upper()  # drop the __ in front of the category
                    # enter to htseq_count_dict
                    library_metadata_dict.setdefault(htseq_count_metadata_category, int(line[1]))
                    # iterate
                line = next(htseq_file_reversed)
        except StopIteration:
            pass

        # error check gene count
        try:
            if crypto_protein_coding_count == 0:
                raise ValueError('NoGeneCountsDetected')
        except ValueError:  #TODO: make this not a static method and add logger
            print('No lines starting with CKF44 have gene counts')

        # rename some key/value pairs
        library_metadata_dict['NOT_ALIGNED_TOTAL'] = library_metadata_dict.pop('NOT_ALIGNED')
        library_metadata_dict['FEATURE_ALIGN_NOT_UNIQUE'] = library_metadata_dict.pop('ALIGNMENT_NOT_UNIQUE')
        library_metadata_dict['AMBIGUOUS_FEATURE'] = library_metadata_dict.pop('AMBIGUOUS')

        # add PROTEIN_CODING_COUNTED
        library_metadata_dict['PROTEIN_CODING_COUNTED'] = crypto_protein_coding_count

        htseq_file.close()
        return library_metadata_dict

    def uniqueAmbiguousProteinCodingCount(self, fastq_simplename):
        """
            intersect bed file with protein coding coords with an alignment file with the htseq annotations added (in this case, grep out __ambiguous)
            :params fastq_simplename: fastq filename minus any path and extention eg /path/to/my_reads_R1.fastq.gz would be my_reads_R1
            :returns: the number of reads (lines) aligning to protein coding coordinates
        """
        #TODO: MAKE THE TRY/CATCH STATEMENT BELOW A FUNCTION AND REPLACE ALL INSTANCES WITH FUNCTION CALL

        # extract bam file from bam_file_list using fastq simple name
        try:
            bam_file = [bam_file for bam_file in self.bam_file_list if fastq_simplename in bam_file][0]
        except IndexError:
            self.logger.debug('bam file not found for %s' % fastq_simplename)  # TODO: improve this logging

        # grep out lines that are ambiguous and map to CKF44 reads (and only CKF44 -- no CNAG)
        extract_ambiguous_protein_coding_reads_cmd = 'samtools view %s| grep ambiguous | grep CKF44| grep -v CNAG| wc -l' %bam_file
        ambiguous_protein_coding_count = int(subprocess.getoutput(extract_ambiguous_protein_coding_reads_cmd))

        return ambiguous_protein_coding_count

    def calculateIntergenicCoverage(self, qual_assess_df):
        """
            calculate coverage of regions in genome between features
            Requires that the number of bases in the intergenic regions of the genome is present in OrganismData_config.ini
            in genome_files/<organism>
        """
        # create new column in qual_assess_df
        qual_assess_df['INTERGENIC_COVERAGE'] = None

        # set up ConfigParser to read OrganismData_config.ini file (this is in each subdir of genome_files)
        kn99_config = configparser.ConfigParser()
        kn99_config.read(os.path.join(self.genome_files, 'KN99', 'OrganismData_config.ini'))
        kn99_config = kn99_config['OrganismData']

        try:
            for index, row in qual_assess_df.iterrows():
                print('...Assessing intergenic coverage for %s' %row['FASTQFILENAME'])
                genotype = list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME']+'.fastq.gz')]['genotype'])[0]
                # set total_intergenic_bases and intergenic_region_bed by organism
                # common error message for setting total_intergenic_bases. Requires %(attribute_name, organism) both strings eg %('kn99_total_intergenic_bases', 'KN99')
                error_msg = 'No attribute %s. Check OrganismData_config.ini in %s'
                if genotype.startswith('CNAG'):
                    try:
                        intergenic_region_bed_path = os.path.join(self.genome_files, 'KN99', kn99_config['intergenic_region_bed'])
                        total_intergenic_bases = int(kn99_config['total_intergenic_bases'])
                    except KeyError:
                        kn99_error_msg = error_msg %('kn99_total_intergenic_bases', 'KN99')
                        self.logger.critical(kn99_error_msg)
                        print(kn99_error_msg)

                # error check intergenic_region_bed_path
                try:
                    if not os.path.isfile(intergenic_region_bed_path):
                        raise FileNotFoundError('IntergenicRegionBedDoesNotExist')
                except FileNotFoundError:
                    intergenic_region_bed_error_msg = 'Intergenic region bed file does not exist at: %s' %intergenic_region_bed_path
                    self.logger.critical(intergenic_region_bed_error_msg)
                    print(intergenic_region_bed_error_msg)

                # extract intergenic_bases_covered from the bam file
                try:
                    bam_file = [bam_file for bam_file in self.bam_file_list if str(row['FASTQFILENAME']) in bam_file][0]
                except IndexError:
                    self.logger.debug('bam file not found for %s' % str(row['FASTQFILENAME']))  # TODO: improve this logging
                intergenic_bases_covered_cmd = 'samtools depth -a -b %s %s | cut -f3 | grep -v 0 | wc -l' %(intergenic_region_bed_path, bam_file)
                num_intergenic_bases_covered = int(subprocess.getoutput(intergenic_bases_covered_cmd))
                qual_assess_df.loc[index, 'INTERGENIC_COVERAGE'] = num_intergenic_bases_covered / float(total_intergenic_bases)

        except NameError:
            self.logger.critical('Cannot calculate INTERGENIC COVERAGE -- total_intergenic_bases or intergenic bed file not found as attribute for organism. Check genome_files/subdirs and each OrganismData_config.ini')
        # return qual_assess_df with INTERGENIC_COVERAGE added
        return qual_assess_df

    def calculateExonicCoverage(self, bam_file, output_directory):
        """
            calculate coverage of exon regions. deposits file in output directory named utils.pathBaseName(bam_file)+'_exonic_coverage.csv'
            :param bam_file: path to bam file
            :param output_directory: path to output directory
        """
        # create df with two columns -- fastqFileName and EXONIC_COVERAGE
        exonic_df = pd.DataFrame()
        exonic_df['fastqFileName'] = [(bam_file)]
        exonic_df['EXONIC_COVERAGE'] = None

        # set up ConfigParser to read OrganismData_config.ini file (this is in each subdir of genome_files)
        kn99_config = configparser.ConfigParser()
        kn99_config.read(os.path.join(self.genome_files, 'KN99', 'OrganismData_config.ini'))
        kn99_config = kn99_config['OrganismData']

        try:
            for index,row in exonic_df.iterrows():
                bam_file = row['fastqFileName']
                print('...Assessing exonic coverage for %s' %bam_file)
                error_msg = 'No attribute %s. Check OrganismData_config.ini in %s'
                # get exon info from config file
                try:
                    exon_region_bed_path = os.path.join(self.genome_files, 'KN99', kn99_config['exon_region_bed'])
                    total_exon_bases = int(kn99_config['total_exon_bases'])
                except KeyError:
                    kn99_error_msg = error_msg %('kn99_total_intergenic_bases', 'KN99')
                    self.logger.critical(kn99_error_msg)
                    print(kn99_error_msg)

                # error check intergenic_region_bed_path
                try:
                    if not os.path.isfile(exon_region_bed_path):
                        raise FileNotFoundError('IntergenicRegionBedDoesNotExist')
                except FileNotFoundError:
                    exonic_region_bed_path_error_msg = 'Intergenic region bed file does not exist at: %s' %exon_region_bed_path
                    self.logger.critical(exonic_region_bed_path_error_msg)
                    print(exonic_region_bed_path_error_msg)

                # extract exonic bases covered by at least one read
                exonic_bases_covered_cmd = 'samtools depth -a -b %s %s | cut -f3 | grep -v 0 | wc -l' %(exon_region_bed_path, bam_file)
                num_exonic_bases_covered = int(subprocess.getoutput(exonic_bases_covered_cmd))

                # add to the df
                exonic_df.loc[index, 'EXONIC_COVERAGE'] = num_exonic_bases_covered / float(total_exon_bases)

        except NameError:
            self.logger.critical('Cannot calculate INTERGENIC COVERAGE -- total_intergenic_bases or intergenic bed file not found as attribute for organism. Check genome_files/subdirs and each OrganismData_config.ini')
        # return qual_assess_df with INTERGENIC_COVERAGE added

        exonic_df['fastqFileName'] = exonic_df['fastqFileName'].apply(lambda x: utils.pathBaseName(x))
        # write
        output_path = os.path.join(output_directory, utils.pathBaseName(bam_file)+'_exonic_coverage.csv')
        exonic_df.to_csv(output_path, index=False)

    def setCryptoGenotypeList(self):
        """ # TODO: Move this to DatabaseObject
            from database_df extract crypto genotype_list (no wildtype CNAG_00000)
        """
        try:
            self.genotype_list = self.standardized_database_df.GENOTYPE.unique()
        except AttributeError:
            print('QualityAssessObject has no standardized_database_df attribute')

    def parseCryptoGenotypeList(self):
        """
            parse a list of genotypes extracted from database_df into ko_gene_list (list of knockout genes) and
            overexpress_gene_list (overexpression genes)
        """
        if len(self.genotype_list) == 0:
            raise IndexError('GenotypeListEmpty')
        for genotype in self.genotype_list:
            if not genotype == 'CNAG_00000':
                if '_over' in genotype:
                    self.overexpress_gene_list.append(genotype.split('_over')[0])
                else:
                    if '.' in genotype:
                        self.ko_gene_list.append(genotype.split('.'))
                    else:
                        self.ko_gene_list.append(genotype)

    def perturbedCheck(self):
        """
           calculate gene coverage of genes in the 'genotype' column of the query_df that do not have suffix _over
           split genotype into two columns, genotype1, genotype2 to address double KO
        """
        # error check if bam_file_list is set
        try:
            if not hasattr(self, 'bam_file_list'):
                raise AttributeError('NoBamFileList')
        except AttributeError:
            print('QualityAssessmentObject does not have attribute bam file list')
        # log which run numbers are present
        self.logger.info('\nThe run numbers in the sheet are: %s' % self.query_df['runNumber'].unique())
        # create genotype_df from query_df[['fastqFileName' and 'genotype']]
        genotype_df = self.query_df[['fastqFileName', 'genotype']]
        # reduce fastqFileName to only simple basename (eg /path/to/some_fastq_R1_001.fastq.gz --> some_fastq_R1_001
        genotype_df['fastqFileName'] = genotype_df['fastqFileName'].apply(lambda x: utils.pathBaseName(x))
        # new column perturbation, if _over in genotype, put 'over' otherwise 'ko'
        genotype_df['pertubation'] = ['over' if '_over' in genotype_column else 'ko' for genotype_column in
                                      genotype_df['genotype']]
        # remove _over from genotype
        genotype_df['genotype'] = genotype_df['genotype'].str.replace('_over', '')
        # split genotype on period. rename column 2 genotype2 if exists. if not, add genotype_2 with values None
        genotype_columns = genotype_df['genotype'].str.split('.', expand=True)
        if len(list(genotype_columns.columns)) == 2:
            genotype_columns.rename(columns={0: 'genotype_1', 1: 'genotype_2'}, inplace=True)
        else:
            genotype_columns.rename(columns={0: 'genotype_1'}, inplace=True)
            genotype_columns['genotype_2'] = None
        # bind genotype_df to genotype_columns
        genotype_df.drop(columns=['genotype'], inplace=True)
        genotype_df = pd.concat([genotype_df, genotype_columns], axis=1)

        # set genomes (this is directly from OrganismData_config.ini in each genome_files subdir)
        s288c_r64_genome = os.path.join(self.genome_files, 'S288C_R64', 'S288C_R64.gff')
        if not os.path.isfile(s288c_r64_genome):
            self.logger.critical('S288C_R64 genome path not valid: %s' % s288c_r64_genome)
            raise FileNotFoundError
        kn99_genome = os.path.join(self.genome_files, 'KN99', 'KN99_stranded_annotations_fungidb_augment.gff')
        if not os.path.isfile(kn99_genome):
            self.logger.critical('S288C_R64 genome path not valid: %s' % kn99_genome)
            raise FileNotFoundError
        h99_genome = os.path.join(self.genome_files, 'H99', 'crNeoH99.gtf')
        if not os.path.isfile(h99_genome):
            self.logger.critical('S288C_R64 genome path not valid: %s' % h99_genome)
            raise FileNotFoundError

        # create columns genotype_1_coverage and genotype_2_coverage
        genotype_df['genotype_1_coverage'] = None
        genotype_df['genotype_2_coverage'] = None
        # loop over rows, calculating coverage for each genotype (testing wither genotype2 is none and perturbation is _over
        for index, row in genotype_df.iterrows():
            if not row['genotype_1'] == 'CNAG_00000':
                # simple name is like this: run_673_s_4_withindex_sequence_TGAGGTT (no containing directories, no extention)
                fastq_simple_name = utils.pathBaseName(row['fastqFileName'])
                genotype_1 = row['genotype_1']
                genotype_2 = row['genotype_2']
                # determine which genome to use -- if CNAG, use kn99
                if genotype_1.startswith('CNAG'):
                    genome = kn99_genome
                    # replace CNAG with CKF44 (in past version of pipeline, KN99 genes were labelled with H99 names. NCBI required change to CKF. Numbering/order is same -- just need to switch CNAG to CKF44)
                    genotype_1 = genotype_1.replace('CNAG', 'CKF44')
                    if genotype_2 is not None and genotype_2.startswith('CNAG'):
                        genotype_2 = genotype_2.replace('CNAG', 'CKF44')
                # for now, if the gene does not start with CNAG, assume it is yeast. TODO: This needs to be changes
                else:
                    genome = s288c_r64_genome
                # get bam files which correspond to query_df fastqFileNames
                try:
                    bam_file = [bam_file for bam_file in self.bam_file_list if fastq_simple_name in bam_file][0]
                except IndexError:
                    self.logger.info('bam file not found for %s' % fastq_simple_name)  # TODO: improve this logging
                    continue
                print('...checking coverage of %s %s' % (genotype_1, genotype_2))
                # extract number of bases in CDS of given gene. Credit: https://www.biostars.org/p/68283/#390427
                num_bases_in_cds_cmd = "grep %s %s | grep CDS | bedtools merge | awk -F\'\t\' \'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}\'" % (
                genotype_1, genome)
                num_bases_in_cds = int(subprocess.getoutput(num_bases_in_cds_cmd))
                # extract number of bases with depth != 0
                num_bases_depth_not_zero_cmd = "grep %s %s | grep CDS | gff2bed | samtools depth -a -b - %s | cut -f3 | grep -v 0 | wc -l" % (
                genotype_1, genome, bam_file)
                num_bases_in_cds_with_one_or_more_read = int(subprocess.getoutput(num_bases_depth_not_zero_cmd))

                # coverage of gene with depth > 1 == (num bases in cds with one or more read) / (number of bases in CDS)
                genotype_df.loc[index, 'genotype_1_coverage'] = num_bases_in_cds_with_one_or_more_read / float(
                    num_bases_in_cds)
                # do the same for genotype_2 if it exists
                if genotype_2 is not None:
                    num_bases_in_cds_cmd = "grep %s %s | grep CDS | bedtools merge | awk -F\'\t\' \'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}\'" % (
                    genotype_2, genome)
                    num_bases_in_cds = int(subprocess.getoutput(num_bases_in_cds_cmd))
                    num_bases_depth_not_zero_cmd = "grep %s %s | grep CDS | gff2bed | samtools depth -a -b - %s | cut -f3 | grep -v 0 | wc -l" % (
                    genotype_2, genome, bam_file)
                    num_bases_in_cds_with_one_or_more_read = int(subprocess.getoutput(num_bases_depth_not_zero_cmd))
                    genotype_df.loc[index, 'genotype_2_coverage'] = num_bases_in_cds_with_one_or_more_read / float(
                        num_bases_in_cds)

        # return genotype check
        genotype_df.columns = [column_name.upper() for column_name in genotype_df.columns]
        return genotype_df[['FASTQFILENAME', 'GENOTYPE_1_COVERAGE', 'GENOTYPE_2_COVERAGE']]

    def ineffectivePerturbationIgvShot(self):
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
