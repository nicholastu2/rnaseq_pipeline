from rnaseq_tools.StandardDataObject import StandardData
from rnaseq_tools.DatabaseObject import DatabaseObject


class QualityAssessmentObject(StandardData):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = ['quality_assessment_filename']
        # TODO: This is a messy and repetitive way of adding expected attributes from children of OrganismData to add to StandardData
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(QualityAssessmentObject, self).__init__(self._add_expected_attributes, **kwargs)
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'QualityAssessmentObject'

    def cryptoPerturbationGenotypeCheck(self):
        # check if necessary paths are entered as attributes
        if not hasattr(self, 'query_sheet_path'):
            raise AttributeError('NoQuerySheetPath')
        if not hasattr(self, 'log2_cpm_path'):
            raise AttributeError('NoLog2CpmPath')
        if not hasattr(self, 'experiment_columns'):
            raise AttributeError('NoExperimentColumns')

        # if align_count_path is present, this is for quality_assess_1 and the prefix should be set as such in column COUNTFILENAME
        if hasattr(self, 'align_count_path'):
            self.standardized_query_df = DatabaseObject.standardizeDatabaseDataframe(self.query_sheet_path,
                                                                                     self.align_count_path)
        # else, no prefix in COUNTFILENAME (see DatabaseObject.StandardizeDatabaseDataframe)
        else:
            self.standardized_query_df = DatabaseObject.standardizeDatabaseDataframe(self.query_sheet_path)

        # TODO: THIS DOES NOT CURRECTLY HANDLE DOUBLE OVER EXPRESSION OR ONE OVER ONE UNDER
        genotype_list = self.standardized_query_df.GENOTYPE.unique()
        ko_gene_list = [] # eg ['CNAG_01020' , ['CNAG_39392','CNAG_48382'], 'CNAG_23421'] where the center item is a double ko
        overexpress_gene_list = [] # expecting no nested lists in this
        for genotype in genotype_list:
            if not genotype == 'CNAG_00000':
                if '_over' in genotype:
                    overexpress_gene_list.append(genotype.split('_over')[0])
                else:
                    if '.' in genotype:
                        ko_gene_list.append(genotype.split('.'))
                    else:
                        ko_gene_list.append(genotype)


    def updateStatusColumn(self):
        pass
