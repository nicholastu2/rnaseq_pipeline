import configparser
from rnaseq_tools import utils
from rnaseq_tools.S288C_R64QualityAssessmentObject import S288C_R64QualityAssessmentObject
import numpy as np

class S288C_R54QualAssessAuditObject(S288C_R64QualityAssessmentObject):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = []
        # This is a method of adding expected attributes to StandardData from StandardData children
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(S288C_R54QualAssessAuditObject, self).__init__(self._add_expected_attributes, **kwargs)
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'CryptoQualAssessAuditObject'
        # create logger
        self.logger = utils.createStandardObjectChildLogger(self, __name__)

        # extract threshold/status from config file #TODO: move to constructor
        qual_assess_config = configparser.ConfigParser()
        qual_assess_config.read(self.config_file)
        qual_assess_1_dict = qual_assess_config['S288C_R64QualityAssessOne']

        # extract thresholds #TODO: CLEAN UP TO DICTIONARY, AUTOMATICALLY EXTRACT
        self.library_size_threshold = int(qual_assess_1_dict['LIBRARY_SIZE_THRESHOLD'])
        self.not_aligned_total_percent_threshold = float(qual_assess_1_dict['NOT_ALIGNED_TOTAL_PERCENT_THRESHOLD'])

        # extract status
        self.library_size_bit_status = int(qual_assess_1_dict['LIBRARY_SIZE_STATUS'])
        self.not_aligned_total_percent_bit_status = int(qual_assess_1_dict['NOT_ALIGNED_TOTAL_PERCENT_STATUS'])

        self.auditQualAssessDataframe()

    def auditQualAssessDataframe(self):
        """
            use rnaseq_pipeline/config/quality_assess_config.ini entries to add status, auto_audit columns
            :params qual_assess_df: a quality_assess_df created by qual_assess_1 or a path to one
            :params bam_file_list: the bam_file_list used in qual_assess_1 -- TODO: make this a list of log2cpm files
            :returns: qual_assess_df with added status and auto_audit columns
        """
        # instantiate a list to hold the status (sum of the status' values when a threshold is failed)
        status_column_list = []

        # loop over the rows of qual_assess_df
        for index, row in self.qual_assess_df.iterrows():
            # extract sample information
            fastq_simple_name = str(row['FASTQFILENAME'])
            # note: only first genotype extracted
            genotype = [self.extractInfoFromQuerySheet(fastq_simple_name, 'genotype1')]

            # extract quality_assessment_metrics
            library_size = int(row['LIBRARY_SIZE'])
            not_aligned_total_percent = float(row['NOT_ALIGNED_TOTAL_PERCENT'])

            # set status_total to 0
            status_total = 0

            # check values against thresholds, add status to flag
            if library_size < self.library_size_threshold:
                status_total += self.library_size_bit_status
            if not_aligned_total_percent > self.not_aligned_total_percent_threshold:
                status_total += self.not_aligned_total_percent_bit_status

            status_column_list.append(status_total)

        self.qual_assess_df['STATUS'] = status_column_list
        self.qual_assess_df['AUTO_AUDIT'] = np.where(self.qual_assess_df.STATUS > 0, 1, 0)
        self.qual_assess_df['STATUS_DECOMP'] = None
        # add status decomposition
        for index, row in self.qual_assess_df.iterrows():
            status = int(row['STATUS'])
            self.qual_assess_df.loc[index, 'STATUS_DECOMP'] = str(utils.decomposeStatus2Bit(status))

        return self.qual_assess_df
