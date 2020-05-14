import os
from rnaseq_tools import utils
from rnaseq_tools.StandardDataObject import StandardData
from rnaseq_tools.DatabaseObject import DatabaseObject
from glob import glob
import pandas as pd


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

        # set optional kwarg arguments # TODO: clean this up
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
            self.ko_gene_list = []  # eg ['CNAG_01020' , ['CNAG_39392','CNAG_48382'], 'CNAG_23421'] where the center item is a double ko
        try:
            self.overexpress_gene_list = kwargs['overexpress_gene_list']
        except KeyError:
            self.overexpress_gene_list = []  # expecting no nested lists in this
        try:
            self.align_count_path = kwargs['align_count_path']
        except KeyError:
            pass

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
            if not hasattr(self, 'query_df'):
                if not os.path.isfile(self.query_path):
                    raise FileNotFoundError('QueryPathNotValid')
                query_df = utils.readInDataframe(self.query_path)
                self.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(query_df)
        except AttributeError:
            print('no standardized query df provided')
        except FileNotFoundError:
            print('query path not valid')
        try:
            if not hasattr(self, 'align_count_path'):
                raise AttributeError('NoAlignCountsPath')
        except AttributeError:
            print('You must pass a path to a directory with alignment files\n'
                  '(typically either the output of align_counts.py or create_experiment.py')

        # create filter (boolean column, used in following line)
        df_wt_filter = self.standardized_database_df['GENOTYPE'] != 'CNAG_00000'
        perturbed_sample_list = list(self.standardized_database_df[df_wt_filter]['COUNTFILENAME'])
        # create path to new sbatch script
        sbatch_job_script_path = os.path.join(self.job_scripts,
                                              'coverage_%s_%s.sbatch' % (self.year_month_day, utils.hourMinuteSecond()))
        # write sbatch script
        print('...writing coverage check sbatch script')
        with open(sbatch_job_script_path, 'w') as sbatch_file:
            sbatch_file.write("#!/bin/bash\n")
            sbatch_file.write("#SBATCH --mem=5G\n")
            sbatch_file.write("#SBATCH -D %s\n" % self.user_rnaseq_pipeline_directory)
            sbatch_file.write("#SBATCH -o sbatch_log/coverage_calculation_%A_%a.out\n")
            sbatch_file.write("#SBATCH -e sbatch_log/coverage_calculation_%A_%a.err\n")
            sbatch_file.write("#SBATCH -J coverage_calculation\n\n")
            sbatch_file.write("ml bedtools\n\n")
            for sample in perturbed_sample_list:
                sorted_alignment_file = sample.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
                sorted_alignment_path = os.path.join(self.align_count_path, sorted_alignment_file)
                coverage_filename = sample.replace('_read_count.tsv', '_coverage.tsv')
                coverage_output_path = os.path.join(self.align_count_path, coverage_filename)
                sbatch_file.write(
                    'bedtools genomecov -ibam %s -bga > %s\n' % (sorted_alignment_path, coverage_output_path))
        print('sbatch script to quantify per base coverage in perturbed samples at %s' % sbatch_job_script_path)
        print('submitting sbatch job. Once this completes, use script quantify_perturbed_coverage.py')
        cmd = 'sbatch %s' % sbatch_job_script_path
        utils.executeSubProcess(cmd)

    def cryptoPertubationBrowserShot(self):
        """
            create IGV browser shot of perturbed gene
        """

    def updateStatusColumn(self):
        pass

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

        try: # TODO: clean this up
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
                return('A sample in the query could not be located in the directory with the alignment files.\n'
                      'Make sure the query provided corresponds to the align_count_path directory')

            output_subdir = os.path.join(self.qorts_output, '[' + bamfilename + ']' + '_qorts')
            utils.mkdirp(output_subdir)
            bamfilename_path = os.path.join(self.align_count_path, bamfilename)
            try:
                if not os.path.exists(bamfilename_path):
                    raise FileExistsError('BamfilePathNotValid')
            except FileExistsError:
                print('path from align_count_path to bamfile does not exist')
            qorts_cmd = 'java -Xmx1G -jar /opt/apps/labs/mblab/software/hartleys-QoRTs-099881f/scripts/QoRTs.jar QC --singleEnded --stranded --keepMultiMapped --generatePlots %s %s %s\n' % (bamfilename_path, self.annotation_file, output_subdir)
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
            qorts_cmd = 'java -Xmx1G -jar /opt/apps/labs/mblab/software/hartleys-QoRTs-099881f/scripts/QoRTs.jar QC --singleEnded --keepMultiMapped --generatePlots %s %s %s\n' % (bamfilename_path, self.annotation_file, output_subdir)
            cmd_list.append(qorts_cmd)

        with open(self.sbatch_script, 'w') as sbatch_file:
            sbatch_file.write('\n')
            sbatch_file.write('\n')

            for cmd in cmd_list:
                sbatch_file.write(cmd)