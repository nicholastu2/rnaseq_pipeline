import os
from rnaseq_tools import utils
from rnaseq_tools.StandardDataObject import StandardData
from rnaseq_tools.DatabaseObject import DatabaseObject
from glob import glob
import pandas as pd
import re

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

    @staticmethod
    def compileData(dir_path, suffix_list):  # TODO: clean up this, parseAlignmentLog and parseCountFile
        """
        get a list of the filenames in the run_#### file that correspond to a given type
        :param dir_path: path to the run_#### directory, generally (and intended to be) in /scratch/mblab/$USER/rnaseq_pipeline/reports
        :param suffix_list: the type of file either novoalign or _read_count (htseq count output) currently set
        :returns: a dataframe containing the files according to their suffix
        """
        # instantiate dataframe
        align_df = pd.DataFrame()
        htseq_count_df = pd.DataFrame()
        # assemble qual_assess_df dataframe
        for suffix in suffix_list: # TODO: error checking: alignment must come before read_count
            # extract files in directory with given suffix
            file_paths = glob("{}/*{}".format(dir_path, suffix))
            for file_path in file_paths:
                # extract fastq filename
                fastq_filename = re.findall(r'(.+?)%s' % suffix, os.path.basename(file_path))[0]
                # set sample name in library_metadata_dict
                library_metadata_dict = {"FASTQFILENAME": fastq_filename}
                if "novoalign" in suffix:
                    library_metadata_dict.update(QualityAssessmentObject.parseAlignmentLog(file_path))
                    align_df = align_df.append(pd.Series(library_metadata_dict), ignore_index=True)
                elif "read_count" in suffix:
                    library_metadata_dict.update(QualityAssessmentObject.parseGeneCount(file_path))
                    htseq_count_df = htseq_count_df.append(pd.Series(library_metadata_dict), ignore_index=True)
        # create dataframe from library_metadata_dict
        qual_assess_df = pd.merge(align_df, htseq_count_df, on='FASTQFILENAME')

        # reformat qual_assess_1 columns
        if "_novoalign.log" in suffix_list and "_read_count.tsv" in suffix_list:  # TODO: CLEAN UP FOLLOWING LINES TO CREATE PERCENT COLUMNS IN SINGLE LINE
            # with_feature_ratio is UNIQUE_ALIGNMENT (from novoalign log) - NO_FEATURE (from htseq-counts output) / no_feature
            # store this for use in htseq column reformatting below
            with_feature_count = (qual_assess_df['UNIQUE_ALIGNMENT'] - qual_assess_df['NO_FEATURE'] - qual_assess_df['AMBIGUOUS_FEATURE'] -
                                  qual_assess_df['TOO_LOW_AQUAL'])
            # create new column TOTAL_ALIGNMENT_PCT as unique + multi maps as percentage of library size
            qual_assess_df['TOTAL_ALIGNMENT'] = (qual_assess_df['UNIQUE_ALIGNMENT'] + qual_assess_df['MULTI_MAP']) / qual_assess_df[
                'LIBRARY_SIZE'].astype('float')
            # total_alignment_count column for use with htseq columns below
            total_alignment_count_column = (qual_assess_df['TOTAL_ALIGNMENT'] * qual_assess_df['LIBRARY_SIZE']).astype('float')
            # convert UNIQUE_ALIGNMENT to percent of library size
            qual_assess_df['UNIQUE_ALIGNMENT'] = qual_assess_df['UNIQUE_ALIGNMENT'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
            # convert MULTI_MAP to percentage of library_size
            qual_assess_df['MULTI_MAP'] = qual_assess_df['MULTI_MAP'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
            # convert NO_MAP tp percentage of library_size
            qual_assess_df['NO_MAP'] = qual_assess_df['NO_MAP'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
            # convert HOMOPOLY_FILTER to percent of library size
            qual_assess_df['HOMOPOLY_FILTER'] = qual_assess_df['HOMOPOLY_FILTER'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
            # convert READ_LENGTH_FILTER to percent of library size
            qual_assess_df['READ_LENGTH_FILTER'] = qual_assess_df['READ_LENGTH_FILTER'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
            # see first line of this if statement for explanation of with_feature_count
            qual_assess_df['WITH_FEATURE_RATIO'] = with_feature_count / qual_assess_df['NO_FEATURE'].astype('float')
            qual_assess_df['WITH_FEATURE'] = with_feature_count / total_alignment_count_column
            qual_assess_df['NO_FEATURE'] = qual_assess_df['NO_FEATURE'] / total_alignment_count_column
            qual_assess_df['FEATURE_ALIGN_NOT_UNIQUE'] = qual_assess_df['FEATURE_ALIGN_NOT_UNIQUE'] / total_alignment_count_column
            qual_assess_df['NOT_ALIGNED_TO_FEATURE'] = qual_assess_df['NOT_ALIGNED_TO_FEATURE'] / total_alignment_count_column
            qual_assess_df['AMBIGUOUS_FEATURE'] = qual_assess_df['AMBIGUOUS_FEATURE'] / total_alignment_count_column
            qual_assess_df['TOO_LOW_AQUAL'] = qual_assess_df['TOO_LOW_AQUAL'] / total_alignment_count_column

        return qual_assess_df.set_index("FASTQFILENAME")

    @staticmethod
    def parseAlignmentLog(alignment_log_file_path):
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
                print('problem with file %s' %alignment_log_file_path)
            # check that the value is both an int and not 0
            if not (isinstance(extracted_value, int) and extracted_value == 0):
                library_metadata_dict.setdefault(alignment_category, extracted_value)
            else:
                print('cannot find %s in %s' %(alignment_category, alignment_log_file_path))

        # close the alignment_file and return
        alignment_file.close()
        return library_metadata_dict

    @staticmethod
    def parseGeneCount(htseq_counts_path):
        """
        count the gene counts that mapped either to genes (see COUNT_VARS at top of script for other features)
        :param htseq_counts_path: a path to a  _read_count.tsv file (htseq-counts output)
        :returns: a dictionary with the keys ALIGNMENT_NOT_UNIQUE, TOO_LOW_AQUAL, AMBIGUOUS_FEATURE, NO_FEATURE
        """
        library_metadata_dict = {}
        # TODO: error checking on keys
        htseq_file = open(htseq_counts_path, 'r')
        htseq_file_reversed = reversed(htseq_file.readlines())

        line = next(htseq_file_reversed)
        while not line.startswith('CNAG'):
            # strip newchar, split on tab
            line = line.strip().split('\t')
            # extract the category of metadata count (eg __alignment_not_unique --> ALIGNMENT_NOT_UNIQUE)
            htseq_count_metadata_category = line[0][2:].upper() # drop the __ in front of the category
            # enter to htseq_count_dict
            library_metadata_dict.setdefault(htseq_count_metadata_category, int(line[1]))
            # iterate
            line = next(htseq_file_reversed)
        # rename some key/value pairs
        library_metadata_dict['NOT_ALIGNED_TO_FEATURE'] = library_metadata_dict.pop('NOT_ALIGNED')
        library_metadata_dict['FEATURE_ALIGN_NOT_UNIQUE'] = library_metadata_dict.pop('ALIGNMENT_NOT_UNIQUE')
        library_metadata_dict['AMBIGUOUS_FEATURE'] = library_metadata_dict.pop('AMBIGUOUS')

        htseq_file.close()
        return library_metadata_dict

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
                sbatch_file.write('\n')
