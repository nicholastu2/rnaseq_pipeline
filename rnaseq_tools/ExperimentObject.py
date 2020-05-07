import pandas as pd
import numpy as np
from rnaseq_tools import utils

class CreateDesignMatrixColumns:
    """
        The intention of this is to create a ExperimentDesign class that is easily extensible to the different varieties of experiments
        we do in the lab. This is a rough (and explicit, as opposed to creating an interface and extending it) draft.
        Currently it is only set up for Zev timecourse experiments
    """

    def __init__(self, experimental_column_headers, contrast_column_header, quality_summary_df, control_value):
        # store cmd line input as class attributes
        self.experimental_column_headers = experimental_column_headers # eg [GENOTYPE:GCN4 & TREATMENT:ESTRADIOL]TIMEPOINT:-1 & 10) GENOTYPE & TREATMENT are the experimental column headers (column names of the query sheet)
        self.contrast_column_header = contrast_column_header # eg in the example above, TIMEPOINT is the contrast heading in the query sheet
        self.quality_summary_df = quality_summary_df
        self.design_df_seed = self.createDesignMatrixSeed()  # remove columns other than experimental + contrast columns from sample_summary_df
        try:
            self.control_value = control_value # ie -1 for TIMEPOINT
            # create experimental controls and contrast groups
            self.experimental_conditions_dict = {}
            self.experimental_conditions_iterable = self.createExperimentalConditionTuples()  # tuples in form ('GCN4', 'Estradiol') in case of [GENOTYPE, TREATMENT]
            self.contrast_conditions = self.createContrastConditions()  # this is without the control condition
            # Using the attributes, create the design matrix
            self.design_df = self.completeDesignMatrix()
        except Exception:
            print('Cannot automatically create your design table. However, the beginning of the table will be depositd in your output.'
                  ' Youll need to finish it by hand.')
            self.design_df = self.design_df_seed


    ### end constructor

    # TODO: include crypto/wt experiments. re-write as general 'experiment' class. extend for various types of experiments.

    def createDesignMatrixSeed(self):
        """
            parse out the experimental_condition_columns + contrast_columns
            of the quality_summary_df and return. This will be the beginning of the design matrix
            Returns:
        """
        columns_of_interest = self.experimental_column_headers + self.contrast_column_header

        if not 'GENOTYPE' in columns_of_interest:
            columns_of_interest.insert(0, 'GENOTYPE')
        if not 'REPLICATE' in columns_of_interest:
            columns_of_interest.insert(1, 'REPLICATE')
        if not 'COUNTFILENAME' in columns_of_interest:
            columns_of_interest.insert(2, 'COUNTFILENAME')

        return self.quality_summary_df[columns_of_interest]

    ### end createDesignMatrixSeed()

    def createExperimentalConditionTuples(self):
        """
            Returns: Tuples created by taking cartesian product of items in the unique values of the experimental_columns in the design_df
        """
        for column_heading in self.experimental_column_headers:
            self.experimental_conditions_dict.setdefault(column_heading, []).extend(
                list(pd.unique(self.design_df_seed[column_heading])))

        return utils.product(*self.experimental_conditions_dict.values())

    ### end createExperimentalConditionsTuples()

    def createContrastConditions(self):
        """
            Returns: List of unique values from the contrast condition column, minus the control_value
        """
        contrast_conditions = list(np.unique(self.design_df_seed[self.contrast_column_header].values))
        contrast_conditions = [str(condition) for condition in contrast_conditions]
        contrast_conditions.remove(self.control_value)

        return contrast_conditions

    ### end createContrastConditions()

    def addDesignMatrixColumn(self, experimental_condition_tuple, contrast_value):
        """
            Creates a column heading in style explained in the FILES section of the rnaseq_pipeline github wiki
            Briefly, the heading will be [experimental conditions]contrast condition (eg [GENETYPE:GCN4 & TREATMENT:ESTRADIOL]TIMEPOINT:-1 & 10)
            Args:
            Returns:
        """
        # append to string while next() for some iterable to create column
        column_heading = ''
        # dictionaries to store column: value to use in identifying the rows in the design_df to add 0s and 1s to
        control_column_identifiers = {}
        contrast_column_identifiers = {}
        for i in range(len(self.experimental_column_headers)):
            # build new column header for the design table
            column_heading = column_heading + ' & ' + self.experimental_column_headers[i] + ':' + str(experimental_condition_tuple[i])
            # column: value to dictionaries
            contrast_column_identifiers.setdefault(self.experimental_column_headers[i], []).append(
                experimental_condition_tuple[i])
            control_column_identifiers.setdefault(self.experimental_column_headers[i], []).append(
                experimental_condition_tuple[i])
        # clean up column heading, add contrast group
        column_heading = '[' + column_heading[3:] + ']' + self.contrast_column_header[0] + ':' + str(
            self.control_value) + ' & ' + str(contrast_value)
        # add contrast column: value to dictionaries
        contrast_column_identifiers.setdefault(self.contrast_column_header[0], []).append(contrast_value)
        control_column_identifiers.setdefault(self.contrast_column_header[0], []).append(self.control_value)
        # create new column in df
        self.design_df_seed[column_heading] = ""
        # fill the column where appropriate with 0s and 1
        self.fillDesignMatrixColumn(column_heading, control_column_identifiers, contrast_column_identifiers)

    ### end addDesignMatrixColumn()

    def fillDesignMatrixColumn(self, column_heading, control_column_identifiers, contrast_column_identifiers):
        """
            Fill the appropriate column with either a zero or 1
            Args:
            Returns:
        """
        # create boolean masks to identify which row is the control and which is the contrast. The dictionaries passed in the arguments are
        # used to filter, which returns 1s (boolean T) where there is a match. If the sum across the rows is equal to the number of keys in
        # the dictionary, then this is a positive match for a row in either the control or contrast group and will be assigned a 0/1 in the loop
        # below
        mask_control = self.design_df_seed.isin(control_column_identifiers).sum(axis=1) == len(
            control_column_identifiers)
        mask_contrast = self.design_df_seed.isin(contrast_column_identifiers).sum(axis=1) == len(
            contrast_column_identifiers)
        # fill row in design_df based on masks above
        for index, row in self.design_df_seed.iterrows():
            if mask_control[index] == True:
                self.design_df_seed.loc[index, column_heading] = '0'
            if mask_contrast[index] == True:
                self.design_df_seed.loc[index, column_heading] = '1'
            else:
                pass

    ### end fillDesignMatrixColumnHeading()

    def completeDesignMatrix(self):
        """
            Append design matrix columns to the design matrix seed
            Args:
            Returns:
        """
        for experimental_condition_tuple in self.experimental_conditions_iterable:
            for contrast_value in self.contrast_conditions:
                self.addDesignMatrixColumn(experimental_condition_tuple, contrast_value)
        # TODO: this needs to be cleaned up -- self.design_df should be modified in place while design_df_seed should be left alone.
        return self.design_df_seed

    ### end completeDesignMatrix()


### end CreateDesignMatrix