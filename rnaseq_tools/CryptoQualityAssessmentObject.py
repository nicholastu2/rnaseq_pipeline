import os
import subprocess
import pandas as pd
import configparser
import sys
from rnaseq_tools import utils
from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject
import numpy as np

# turn off SettingWithCopyWarning in pandas
pd.options.mode.chained_assignment = None
#TODO: CURRENTLY, ALWAYS CALCULATES COVERAGES -- CREATE FLAG TO SKIP THIS STEP IF PASSED
#TODO: CALCULATE COVERAGE ONCE, STORE AS BED FILE, USE BED RATHER THAN QUANTIFYING BAM EVERYTIME
class CryptoQualityAssessmentObject(QualityAssessmentObject):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = []
        # This is a method of adding expected attributes to StandardData from StandardData children
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(CryptoQualityAssessmentObject, self).__init__(self._add_expected_attributes, **kwargs)
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'CryptoQualityAssessmentObject'

        # for ordering columns below. genotype_1_coverage and genotype_2_coverage added if coverage_check is passed
        self.column_order = ['FASTQFILENAME','LIBRARY_SIZE', 'EFFECTIVE_LIBRARY_SIZE', 'EFFECTIVE_UNIQUE_ALIGNMENT', 'EFFECTIVE_UNIQUE_ALIGNMENT_PERCENT',
                             'MULTI_MAP_PERCENT', 'PROTEIN_CODING_TOTAL', 'PROTEIN_CODING_TOTAL_PERCENT', 'PROTEIN_CODING_COUNTED',
                             'PROTEIN_CODING_COUNTED_PERCENT', 'AMBIGUOUS_FEATURE_PERCENT', 'NO_FEATURE_PERCENT',
                             'INTERGENIC_COVERAGE', 'NOT_ALIGNED_TOTAL_PERCENT', 'GENOTYPE_1_COVERAGE', 'GENOTYPE_2_COVERAGE',
                             'NAT_COVERAGE', 'G418_COVERAGE', 'OVEREXPRESSION_FOW', 'NO_MAP_PERCENT', 'HOMOPOLY_FILTER_PERCENT', 'READ_LENGTH_FILTER_PERCENT',
                             'TOO_LOW_AQUAL_PERCENT', 'rRNA_PERCENT', 'nctrRNA_PERCENT']
        # set organismdata_config.ini in genome_files
        crypto_genome_files_config_path = os.path.join(self.genome_files, 'KN99', 'OrganismData_config.ini')
        utils.configure(self, crypto_genome_files_config_path, 'OrganismData', os.path.join(self.genome_files, 'KN99'))


    def formatQualAssessDataFrame(self, qual_assess_df):
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

        return qual_assess_df[self.column_order]


    def auditQualAssessDataFrame(self, query_df_path, qual_assess_df, bam_file_list):
        """
            use rnaseq_pipeline/config/quality_assess_config.ini entries to add status, auto_audit columns
            :params qual_assess_df: a quality_assess_df created by qual_assess_1 or a path to one
            :params bam_file_list: the bam_file_list used in qual_assess_1 -- TODO: make this a list of log2cpm files
            :returns: qual_assess_df with added status and auto_audit columns
        """
        # read in median_wt_expression_by_timepoint_treatment_df
        try:
            median_wt_expression_by_timepoint_treatment_df = utils.readInDataframe(self.median_wt_expression_by_timepoint_treatment)
            # SET INDEX ON (gene_id, TREATMENT, TIMEPOINT) note: timepoint is read in as an int
            median_wt_expression_by_timepoint_treatment_df = median_wt_expression_by_timepoint_treatment_df.set_index(['gene_id', 'TREATMENT', 'TIMEPOINT'])
        except AttributeError:
            self.logger.critical('genome files config in constructor did not work')
            print('genome files config in constructor did not work')

        # read in qual_assess_df if necessary
        if os.path.isfile(qual_assess_df):
            qual_assess_df = utils.readInDataframe(qual_assess_df)
        # read in query_df
        if os.path.isfile(query_df_path):
            self.query_df = utils.readInDataframe(query_df_path)
        # remove na/nan fastqFileName rows
        self.query_df = self.query_df[~self.query_df.fastqFileName.isna()]

        # extract threshold/status from config file #TODO: move to constructor
        qual_assess_config = configparser.ConfigParser()
        qual_assess_config.read(self.config_file)
        qual_assess_1_dict = qual_assess_config['KN99QualityAssessOne']

        # extract thresholds
        protein_coding_total_threshold = int(qual_assess_1_dict['PROTEIN_CODING_TOTAL_THRESHOLD'])
        not_aligned_total_percent_threshold = float(qual_assess_1_dict['NOT_ALIGNED_TOTAL_PERCENT_THRESHOLD'])
        perturbed_coverage_threshold = float(qual_assess_1_dict['PERTURBED_COVERAGE_THRESHOLD'])
        expected_nat_coverage_threshold = float(qual_assess_1_dict['EXPECTED_NAT_COVERAGE_THRESHOLD'])
        expected_nat_log2cpm_threshold = float(qual_assess_1_dict['EXPECTED_NAT_LOG2CPM_THRESHOLD'])
        unexpected_nat_coverage_threshold = float(qual_assess_1_dict['UNEXPECTED_NAT_COVERAGE_THRESHOLD'])
        unexpected_nat_log2cpm_threshold = float(qual_assess_1_dict['UNEXPECTED_NAT_LOG2CPM_THRESHOLD'])
        g418_log2cpm_threshold = float(qual_assess_1_dict['G418_LOG2CPM_THRESHOLD'])
        overexpression_fow_threshold = float(qual_assess_1_dict['OVEREXPRESSION_FOW_THRESHOLD'])


        # extract status
        protein_coding_total_bit_status = int(qual_assess_1_dict['PROTEIN_CODING_TOTAL_STATUS'])
        not_aligned_total_percent_bit_status = int(qual_assess_1_dict['NOT_ALIGNED_TOTAL_PERCENT_STATUS'])
        perturbed_coverage_bit_status = int(qual_assess_1_dict['PERTURBED_COVERAGE_STATUS'])
        nat_expected_marker_status = int(qual_assess_1_dict['NAT_EXPECTED_MARKER_STATUS'])
        nat_unexpected_marker_status = int(qual_assess_1_dict['NAT_UNEXPECTED_MARKER_STATUS'])
        g418_expected_marker_status = int(qual_assess_1_dict['G418_EXPECTED_MARKER_STATUS'])
        g418_unexpected_marker_status = int(qual_assess_1_dict['G418_UNEXPECTED_MARKER_STATUS'])
        overexpression_fow_status = int(qual_assess_1_dict['OVEREXPRESSION_FOW_STATUS'])
        no_metadata_marker_status = int(qual_assess_1_dict['NO_METADATA_MARKER_STATUS'])

        # instantiate a list to hold the status (sum of the status' values when a threshold is failed)
        status_column_list = []
        # loop over the rows of qual_assess_df
        for index, row in qual_assess_df.iterrows():
            # extract fastq simple name
            fastq_simple_name = str(row['FASTQFILENAME'])

            # extract treatment # TODO: MAKE EXTRACTING INFO FROM QUERY_DF FROM QUAL_ASSES A FUNCTION IN QUALITYASSESSMENTOBJECT
            try:
                sample_treatment = list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME'] + '.fastq.gz')]['treatment'])[0]
            except AttributeError:
                sys.exit('You must pass a query df')

            # extract timePoint
            sample_timepoint = list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME'] + '.fastq.gz')]['timePoint'])[0]

            # extract genotype
            genotype = list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME'] + '.fastq.gz')]['genotype'])[0]

            # split on period to separate double KO. note that this is now a list, even if one item
            genotype = genotype.split('.')

            # get markers
            marker_1 = str(list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME'] + '.fastq.gz')]['marker_1'])[0])
            marker_2 = str(list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME'] + '.fastq.gz')]['marker_2'])[0])

            # log2cpm path TODO: THIS NEEDS TO BE FIXED -- INPUT LIST OF LOG2_CPM AND REDO NAMING CONVENTION IN LOG2 CPM SCRIPT TO INCLUDE CONTAINING FOLDER
            bam_path = [bam_file for bam_file in bam_file_list if fastq_simple_name in bam_file][0]
            run_dirpath = utils.dirPath(utils.dirPath(bam_path))
            log2_cpm_path = os.path.join(run_dirpath, 'count', 'log2_cpm.csv')
            try:
                if not os.path.isfile(log2_cpm_path):
                    raise FileNotFoundError('Log2CpmSheetNotFound')
            except FileNotFoundError:
                error_msg = 'log2_cpm sheet path %s not valid' %log2_cpm_path
                self.logger.critical(error_msg)
                print(error_msg)

            # set overexpression flag
            overexpression_flag = False
            if '_over' in genotype[0]:
                overexpression_flag = True
                # set perturbed_gene variable
                perturbed_gene = genotype[0].replace('_over', '')
                # extract overexpression log2_cpm
                overexpression_log2cpm = float(self.extractLog2cpm(perturbed_gene, fastq_simple_name, log2_cpm_path))
                # get wildtype log2_cpm from median_wt-expression_by_timepoint_treatment
                try:
                    wt_log2cpm = float(median_wt_expression_by_timepoint_treatment_df.loc[(perturbed_gene, sample_treatment, int(sample_timepoint)), 'MEDIAN_LOG2CPM'])
                except KeyError:
                    perturbed_gene = perturbed_gene.replace('CNAG', 'CKF44')
                    wt_log2cpm = float(median_wt_expression_by_timepoint_treatment_df.loc[
                                           (perturbed_gene, sample_treatment, int(sample_timepoint)), 'MEDIAN_LOG2CPM'])

            # extract quality_assessment_metrics
            library_protein_coding_total = int(row['PROTEIN_CODING_TOTAL'])
            not_aligned_total_percent = float(row['NOT_ALIGNED_TOTAL_PERCENT'])
            if row['GENOTYPE_1_COVERAGE'] is not None:
                library_genotype_1_coverage = float(row['GENOTYPE_1_COVERAGE'])
            else:
                library_genotype_1_coverage = -1
            if row['GENOTYPE_2_COVERAGE'] is not None:
                library_genotype_2_coverage = float(row['GENOTYPE_2_COVERAGE'])
            else:
                library_genotype_2_coverage = -1

            # extract NAT log2cpm
            nat_log2cpm = self.extractLog2cpm('CNAG_NAT', fastq_simple_name, log2_cpm_path)
            # extract NAT coverage
            try:
                nat_coverage = float(row['NAT_COVERAGE'])
            except KeyError:
                self.logger.critical('%s no nat coverage' %row)
            # extract G418 log2cpm
            g418_log2cpm = self.extractLog2cpm('CNAG_G418', fastq_simple_name, log2_cpm_path)

            # set status_total to 0
            status_total = 0

            # check values against thresholds, add status to flag
            if library_protein_coding_total < protein_coding_total_threshold:
                status_total += protein_coding_total_bit_status
            if not_aligned_total_percent > not_aligned_total_percent_threshold:
                status_total += not_aligned_total_percent_bit_status
            # if the overexpression_flag is set, eval based on expression
            if overexpression_flag:
                overexpression_fow = overexpression_log2cpm - wt_log2cpm
                qual_assess_df.loc[index, 'OVEREXPRESSION_FOW'] = overexpression_fow
                if overexpression_fow < overexpression_fow_threshold:
                    status_total += overexpression_fow_status
            # else, evaluate KO based on coverage
            else:
                if library_genotype_1_coverage > perturbed_coverage_threshold or library_genotype_2_coverage > perturbed_coverage_threshold:
                    status_total += perturbed_coverage_bit_status

            # test wildtypes for marker coverage and expression
            if genotype[0] == 'CNAG_00000':
                if nat_coverage > unexpected_nat_coverage_threshold or nat_log2cpm > unexpected_nat_log2cpm_threshold:
                    status_total += nat_unexpected_marker_status
                if g418_log2cpm > g418_log2cpm_threshold:
                    status_total += g418_unexpected_marker_status
            # if a perturbed sample, first check that marker information is present and flag it if it is not
            else:
                if marker_1 == 'nan' or marker_1 is None or (len(genotype) > 1 and marker_2 == 'nan' or marker_2 == 'none'):
                    status_total+=no_metadata_marker_status
                # if perturbed, and marker information is present, test the markers
                else:
                    if marker_1 == 'NAT':
                        if nat_coverage < expected_nat_coverage_threshold or nat_log2cpm < expected_nat_log2cpm_threshold:
                            status_total += nat_expected_marker_status
                        if g418_log2cpm > g418_log2cpm_threshold:
                            status_total += g418_unexpected_marker_status
                    elif marker_1 == 'G418':
                        if g418_log2cpm < g418_log2cpm_threshold:
                            status_total += g418_expected_marker_status
                        if nat_coverage > unexpected_nat_coverage_threshold or nat_log2cpm > unexpected_nat_log2cpm_threshold:
                            status_total += nat_unexpected_marker_status
                    # TODO: THIS NEEDS TO BE FIXED FOR THE INSTANCE IN WHICH A SINGLE MARKER IS NOTED WITH BOTH MARKERS B/C OF MANUAL ENTRY (SEE STRAINS 2274 AND 1351)
                    # TEST THIS WITH THE STRAINS ABOVE
                    if len(genotype) > 1 or marker_2 != 'nan':  # note: unentered 2nd markers for double KO should be caught in the if statement above
                        if marker_2 == 'NAT':
                            if marker_1 == 'NAT':
                                self.logger.critical('%s has two NAT markers in the metadata' %fastq_simple_name)
                            if nat_coverage < expected_nat_coverage_threshold or nat_log2cpm < expected_nat_log2cpm_threshold:
                                status_total += nat_expected_marker_status
                            if g418_log2cpm > g418_log2cpm_threshold:
                                status_total += g418_unexpected_marker_status
                        elif marker_2 == 'G418':
                            if marker_1 == 'G418':
                                self.logger.critical('%s has two G418 markers in the metadata' %fastq_simple_name)
                            if g418_log2cpm < g418_log2cpm_threshold:
                                status_total += g418_expected_marker_status
                            if nat_coverage > unexpected_nat_coverage_threshold or nat_log2cpm > unexpected_nat_log2cpm_threshold:
                                status_total += nat_unexpected_marker_status

            status_column_list.append(status_total)

        qual_assess_df['STATUS'] = status_column_list
        qual_assess_df['AUTO_AUDIT'] = np.where(qual_assess_df.STATUS > 0, 1, 0)
        qual_assess_df['STATUS_DECOMP'] = None
        # add status decomposition
        for index,row in qual_assess_df.iterrows():
            status = int(row['STATUS'])
            qual_assess_df.loc[index, 'STATUS_DECOMP'] = str(utils.decomposeStatus2Bit(status))

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
        genotype_df['perturbation'] = ['over' if '_over' in genotype_column else 'ko' for genotype_column in
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

        # read in KN99 OrganismData_config.ini
        kn99_config_filepath = os.path.join(self.genome_files, 'KN99', 'OrganismData_config.ini')
        if not os.path.isfile(kn99_config_filepath):
            self.logger.critical('kn99 config file not found at %s' %kn99_config_filepath)
            raise FileNotFoundError('KN99ConfigFilpathNotFound')
        config = configparser.ConfigParser()
        config.read(kn99_config_filepath)
        kn99_organism_data_dict = config['OrganismData']
        # extract genome path
        kn99_annotation_path = os.path.join(self.genome_files, 'KN99', kn99_organism_data_dict['annotation_file'])
        if not os.path.isfile(kn99_config_filepath):
            self.logger.critical('kn99 annotation file not found at %s' %kn99_config_filepath)
            raise FileNotFoundError('KN99AnnotationFileNotFound')
        nat_bases_in_cds = int(kn99_organism_data_dict['NAT_cds_length'])
        g418_bases_in_cds = int(kn99_organism_data_dict['G418_cds_length'])

        # set feature over which to take percentage of reads (CDS in this case)
        feature = 'CDS'

        # create columns genotype_1_coverage, genotype_2_coverage, overexpression_fow (fold over wildtype), NAT_coverage, G418_coverage
        genotype_df['genotype_1_coverage'] = None
        genotype_df['genotype_2_coverage'] = None
        genotype_df['overexpression_fow'] = None
        genotype_df['NAT_coverage'] = None
        genotype_df['G418_coverage'] = None

        # loop over rows, calculating coverage for each genotype (testing wither genotype2 is none and perturbation is _over
        for index, row in genotype_df.iterrows():
            # simple name is like this: run_673_s_4_withindex_sequence_TGAGGTT (no containing directories, no extention)
            fastq_simple_name = utils.pathBaseName(row['fastqFileName'])
            # get bam files which correspond to query_df fastqFileNames
            try:
                bam_file = [bam_file for bam_file in self.bam_file_list if fastq_simple_name in bam_file][0]
            except IndexError:
                self.logger.info('bam file not found for %s' % fastq_simple_name)  # TODO: improve this logging
                continue

            # calculate marker coverages
            print('...calculation NAT coverage for %s' %fastq_simple_name)
            genotype_df.loc[index, 'NAT_coverage'] = self.calculatePercentFeatureCoverage(feature, 'CNAG_NAT', kn99_annotation_path, bam_file, nat_bases_in_cds)
            print('...calculation G418 coverage for %s' %fastq_simple_name)
            genotype_df.loc[index, 'G418_coverage'] = self.calculatePercentFeatureCoverage(feature, 'CNAG_G418', kn99_annotation_path, bam_file, g418_bases_in_cds)

            # if perturbed, calculate perturbed gene coverage
            if not row['genotype_1'] == 'CNAG_00000' and not row['perturbation'] == 'over':
                genotype_1 = row['genotype_1']
                genotype_2 = row['genotype_2']
                # determine which genome to use -- if CNAG, use kn99
                if not genotype_1.startswith('CNAG'):
                    raise ValueError('%sNotRecognizedCryptoGenotype' %genotype_1)
                # replace CNAG with CKF44 (in past version of pipeline, KN99 genes were labelled with H99 names. NCBI required change to CKF. Numbering/order is same -- just need to switch CNAG to CKF44)
                genotype_1 = genotype_1.replace('CNAG', 'CKF44')
                if genotype_2 is not None and genotype_2.startswith('CNAG'):
                    genotype_2 = genotype_2.replace('CNAG', 'CKF44')
                print('...checking coverage of %s %s in %s' % (genotype_1, genotype_2, fastq_simple_name))
                genotype_df.loc[index, 'genotype_1_coverage'] = self.calculatePercentFeatureCoverage(feature, genotype_1, kn99_annotation_path, bam_file)
                # do the same for genotype_2 if it exists
                if genotype_2 is not None:
                    genotype_df.loc[index, 'genotype_2_coverage'] = self.calculatePercentFeatureCoverage(feature, genotype_2, kn99_annotation_path, bam_file)

        # return genotype check
        genotype_df.columns = [column_name.upper() for column_name in genotype_df.columns]
        return genotype_df[['FASTQFILENAME', 'GENOTYPE_1_COVERAGE', 'GENOTYPE_2_COVERAGE', 'OVEREXPRESSION_FOW', 'NAT_COVERAGE', 'G418_COVERAGE']]