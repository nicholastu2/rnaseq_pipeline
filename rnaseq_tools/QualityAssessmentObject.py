import sys
import os
from rnaseq_tools import utils
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
        # set standardDirectory structure
        self.standardDirectoryStructure()
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'QualityAssessmentObject'

        # set optional kwarg arguments
        try:
            self.standardized_database_df = kwargs['standardized_database_df']
        except KeyError:
            pass
        try:
            self.genotype_list = kwargs['genotype_list']
        except KeyError:
            self.genotype_list = []
        try:
            self.ko_gene_list = kwargs['ko_gene_list']
        except KeyError:
            self.ko_gene_list = [] # eg ['CNAG_01020' , ['CNAG_39392','CNAG_48382'], 'CNAG_23421'] where the center item is a double ko
        try:
            self.overexpress_gene_list = kwargs['overexpress_gene_list']
        except KeyError:
            self.overexpress_gene_list = [] # expecting no nested lists in this

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


    def cryptoPerturbationExpressionCheck(self):
        """
            check expression against wt expression as mean_log2_perturbed_treatment_timepoint / mean_wt_treatment_timepoint
        """
        # check if necessary paths are entered as attributes
        if not hasattr(self, 'query_sheet_path'):
            raise AttributeError('NoQuerySheetPath')
        if not hasattr(self, 'log2_cpm_path'):
            raise AttributeError('NoLog2CpmPath')
        if not hasattr(self, 'experiment_columns'):
            raise AttributeError('NoExperimentColumns')

        # if standardized_database_df is not yet set, do so
        # if align_count_path is present, this is for quality_assess_1 and the prefix should be set as such in column COUNTFILENAME
        if hasattr(self, 'align_count_path') and not hasattr(self, 'standardized_database_df'):
            self.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(self.query_sheet_path,
                                                                                     self.align_count_path)
        # else, no prefix in COUNTFILENAME (see DatabaseObject.StandardizeDatabaseDataframe)
        elif not hasattr(self, 'standardized_database_df'):
            self.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(self.query_sheet_path)

        for gene in self.ko_gene_list:
            pass
            # extract mean perturbed log2cpm expression and treatment_timepoint
            # if gene expression is greater than 20% of wt in treatment_timempoint, flag, create diagnostic dataframe and browser shot

        for gene in self.overexpress_gene_list:
            pass
            # extract mean perturbed log2cpm expression and treatment_timepoint
            # if gene expression is less than 99% of wt in same treatment_timempoint, flag, create diagnostic dataframe and browser shot

    def coverageCheck(self):
        """
            check coverage of genebody
        """
        print('...extracting perturbed samples from query sheet for coverage check sbatch script')
        # check if all necessary components are present
        try:
            if not hasattr(self, 'query_path'):
                raise AttributeError('NoQueryPath')
            if not hasattr(self, 'query_df')
                query_df = utils.readInDataframe(self.query_path)
            if not hasattr(self, 'standardized_database_df'):
                self.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(query_df)
        except AttributeError:
            sys.exit('no standardized query df provided')

        # extract perturbed samples' COUNTFILENAME
        self.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(qa.query_df)
        # create filter (boolean column, used in following line)
        df_wt_filter = self.standardized_database_df['GENOTYPE'] != 'CNAG_00000'
        perturbed_sample_list = list(self.standardized_database_df[df_wt_filter]['COUNTFILENAME'])
        # create path to new sbatch script
        sbatch_job_script_path = os.path.join(self.job_scripts, 'coverage_%s_%s.sbatch' %(self.year_month_day, utils.hourMinuteSecond()))
        # write sbatch script
        print('...writing coverage check sbatch script')
        with open(sbatch_job_script_path, 'w') as sbatch_file:
            sbatch_file.write("#!/bin/bash\n")
            sbatch_file.write("#SBATCH --mem=5G\n")
            sbatch_file.write("#SBATCH -D %s\n" %self.user_rnaseq_pipeline_directory)
            sbatch_file.write("#SBATCH -o sbatch_log/coverage_calculation_%A_%a.out\n")
            sbatch_file.write("#SBATCH -e sbatch_log/coverage_calculation_%A_%a.err\n")
            sbatch_file.write("#SBATCH -J coverage_calculation\n\n")
            sbatch_file.write("ml bedtools\n\n")
            for sample in perturbed_sample_list:
                sorted_alignment_file = sample.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
                sorted_alignment_path = os.path.join(self.align_count_path, sorted_alignment_file)
                coverage_filename = sample.replace('_read_count.tsv', '_coverage.tsv')
                coverage_output_path = os.path.join(self.align_count_path, coverage_filename)
                sbatch_file.write('bedtools genomecov -ibam %s -bga > %s\n' %(sorted_alignment_path, coverage_output_path))
        print('sbatch script to quantify per base coverage in perturbed samples at %s' %sbatch_job_script_path)
        print('submitting sbatch job. Once this completes, use script quantify_perturbed_coverage.py')
        cmd = 'sbatch %s'%sbatch_job_script_path
        utils.executeSubProcess(cmd)

    def cryptoPertubationBrowserShot(self):
        """
            create IGV browser shot of perturbed gene
        """



    def updateStatusColumn(self):
        pass
