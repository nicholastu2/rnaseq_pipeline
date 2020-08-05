import os
import pandas as pd
from rnaseq_tools import utils
from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject

# turn off SettingWithCopyWarning in pandas
pd.options.mode.chained_assignment = None


# TODO: CURRENTLY, ALWAYS CALCULATES COVERAGES -- CREATE FLAG TO SKIP THIS STEP IF PASSED
# TODO: CALCULATE COVERAGE ONCE, STORE AS BED FILE, USE BED RATHER THAN QUANTIFYING BAM EVERYTIME
class S288C_R64QualityAssessmentObject(QualityAssessmentObject):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = []
        # This is a method of adding expected attributes to StandardData from StandardData children
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(S288C_R64QualityAssessmentObject, self).__init__(self._add_expected_attributes, **kwargs)
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'CryptoQualityAssessmentObject'

        # for ordering columns below. genotype_1_coverage and genotype_2_coverage added if coverage_check is passed
        self.column_order = ['FASTQFILENAME', 'LIBRARY_SIZE', 'EFFECTIVE_LIBRARY_SIZE', 'EFFECTIVE_UNIQUE_ALIGNMENT',
                             'EFFECTIVE_UNIQUE_ALIGNMENT_PERCENT',
                             'MULTI_MAP_PERCENT', 'PROTEIN_CODING_TOTAL', 'PROTEIN_CODING_TOTAL_PERCENT',
                             'PROTEIN_CODING_COUNTED',
                             'PROTEIN_CODING_COUNTED_PERCENT', 'AMBIGUOUS_FEATURE_PERCENT', 'NO_FEATURE_PERCENT',
                             'INTERGENIC_COVERAGE', 'NOT_ALIGNED_TOTAL_PERCENT', 'GENOTYPE_1_COVERAGE',
                             'GENOTYPE_2_COVERAGE',
                             'NAT_COVERAGE', 'G418_COVERAGE', 'OVEREXPRESSION_FOW', 'NO_MAP_PERCENT',
                             'HOMOPOLY_FILTER_PERCENT', 'READ_LENGTH_FILTER_PERCENT',
                             'TOO_LOW_AQUAL_PERCENT', 'rRNA_PERCENT', 'nctrRNA_PERCENT']
        # set organismdata_config.ini in genome_files
        S288C_R642_genome_files_config_path = os.path.join(self.genome_files, 'S288C_R642', 'OrganismData_config.ini')
        utils.configure(self, S288C_R642_genome_files_config_path, 'OrganismData', os.path.join(self.genome_files, 'S288C_R642'))

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
            htseq_count_df = htseq_count_df.append(pd.Series(library_metadata_dict), ignore_index=True)
        print('\nDone parsing count files\n')

        # concat df_list dataframes together on the common column. CREDIT: https://stackoverflow.com/a/56324303/9708266
        qual_assess_df = pd.merge(align_df, htseq_count_df, on='FASTQFILENAME')
        qual_assess_df = self.formatQualAssessDataFrame(qual_assess_df)

        return qual_assess_df

    def parseGeneCount(self, htseq_counts_path):
        """
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
                if line.startswith('--'):
                    # strip newchar, split on tab
                    line = line.strip().split('\t')
                    # extract the category of metadata count (eg __alignment_not_unique --> ALIGNMENT_NOT_UNIQUE)
                    htseq_count_metadata_category = line[0][2:].upper()  # drop the __ in front of the category
                    # enter to htseq_count_dict
                    library_metadata_dict.setdefault(htseq_count_metadata_category, int(line[1]))
                    # iterate
                line = next(htseq_file_reversed)
        except StopIteration:
            pass

        # rename some key/value pairs
        library_metadata_dict['NOT_ALIGNED_TOTAL'] = library_metadata_dict.pop('NOT_ALIGNED')
        library_metadata_dict['FEATURE_ALIGN_NOT_UNIQUE'] = library_metadata_dict.pop('ALIGNMENT_NOT_UNIQUE')
        library_metadata_dict['AMBIGUOUS_FEATURE'] = library_metadata_dict.pop('AMBIGUOUS')

        htseq_file.close()
        return library_metadata_dict
