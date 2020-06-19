import os
import subprocess
import re
import pandas as pd
from glob import glob
from rnaseq_tools import utils
from rnaseq_tools.StandardDataObject import StandardData
from rnaseq_tools.DatabaseObject import DatabaseObject

# turn off SettingWithCopyWarning
pd.options.mode.chained_assignment = None

class QualityAssessmentObject(StandardData):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = ['quality_assessment_filename', 'nextflow_list_of_files', 'coverage_check_flag']
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

        # set optional kwarg arguments # TODO: clean this up
        try:
            self.query_path = kwargs['query_path']
            try:
                self.query_df = utils.readInDataframe(self.query_path)
            except ValueError:
                self.logger.critical('query_path not valid')
            except FileNotFoundError:
                self.logger.critical('%s  --> query_path not valid' %self.query_path)
        except KeyError:
            pass
        try:
            self.nextflow_list_of_files = kwargs['nextflow_list_of_files']
        except KeyError:
            pass
        try:
            self.coverage_check_flag = kwargs['coverage_check_flag']
        except KeyError:
            pass
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
            self.quality_assess_dir_path = kwargs['quality_assess_dir_path']
        except KeyError:
            pass
        # the suffixes of the log/metadata info in list form. eg ['_novoalign.log', '_read_count.tsv'] for novoalign htseq count respectively
        try:
            self.log_suffix_list = kwargs['log_suffix_list']
        except KeyError:
            pass

    def compileData(self):  # TODO: clean up this, parseAlignmentLog and parseCountFile
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
        for suffix in self.log_suffix_list:  # TODO: error checking: alignment must come before read_count
            try:  # TODO: clean this up. ugly method of testing whether a list of filepaths has been passed or not. this is for the nextflow pipeline
                file_paths = self.nextflow_list_of_files
            except AttributeError:
                # extract files in directory with given suffix
                file_paths = glob("{}/*{}".format(self.quality_assess_dir_path, suffix))
            for file_path in file_paths:
                # extract fastq filename
                fastq_basename = re.findall(r'(.+?)%s' % suffix, os.path.basename(file_path))[0]
                # set sample name in library_metadata_dict
                library_metadata_dict = {"FASTQFILENAME": fastq_basename}
                if "novoalign" in suffix:
                    library_metadata_dict.update(QualityAssessmentObject.parseAlignmentLog(file_path))
                    align_df = align_df.append(pd.Series(library_metadata_dict), ignore_index=True)
                elif "read_count" in suffix:
                    library_metadata_dict.update(QualityAssessmentObject.parseGeneCount(file_path))
                    htseq_count_df = htseq_count_df.append(pd.Series(library_metadata_dict), ignore_index=True)
        # concat df_list dataframes together on the common column. CREDIT: https://stackoverflow.com/a/56324303/9708266
        qual_assess_df = pd.merge(align_df, htseq_count_df, on='FASTQFILENAME')

        # reformat qual_assess_1 columns
        if "_novoalign.log" in self.log_suffix_list and "_read_count.tsv" in self.log_suffix_list:  # TODO: CLEAN UP FOLLOWING LINES TO CREATE PERCENT COLUMNS IN SINGLE LINE
            # store library size column
            library_size_column = qual_assess_df['LIBRARY_SIZE'].astype('float')
            # total_with_feature is the unique alignments MINUS N0_FEATURE, AMBIGUOUS_FEATURE and TOO_LOW_AQUAL
            with_feature_column = (qual_assess_df['UNIQUE_ALIGNMENT'] - qual_assess_df['NO_FEATURE'] - qual_assess_df[
                'AMBIGUOUS_FEATURE'] - qual_assess_df['TOO_LOW_AQUAL'])
            unique_alignment_total = qual_assess_df['UNIQUE_ALIGNMENT'].astype('float')

            # create new column TOTAL_ALIGNMENT_PCT as unique + multi maps as percentage of library size
            qual_assess_df['TOTAL_ALIGNMENT'] = (qual_assess_df['UNIQUE_ALIGNMENT'] + qual_assess_df[
                'MULTI_MAP']) / library_size_column
            # convert novoalign columns to percents of library_size_column
            total_alignment_count_column = qual_assess_df['TOTAL_ALIGNMENT'] * library_size_column
            qual_assess_df['UNIQUE_ALIGNMENT'] = qual_assess_df['UNIQUE_ALIGNMENT'] / library_size_column
            qual_assess_df['MULTI_MAP'] = qual_assess_df['MULTI_MAP'] / library_size_column
            qual_assess_df['NO_MAP'] = qual_assess_df['NO_MAP'] / library_size_column
            qual_assess_df['HOMOPOLY_FILTER'] = qual_assess_df['HOMOPOLY_FILTER'] / library_size_column
            qual_assess_df['READ_LENGTH_FILTER'] = qual_assess_df['READ_LENGTH_FILTER'] / library_size_column

            # convert htseq count columns to percent of aligned reads
            qual_assess_df['NOT_ALIGNED_TOTAL'] = qual_assess_df['NOT_ALIGNED_TOTAL'] / library_size_column  # TODO: CHECK THIS
            qual_assess_df['WITH_FEATURE'] = with_feature_column / unique_alignment_total
            qual_assess_df['NO_FEATURE'] = qual_assess_df['NO_FEATURE'] / unique_alignment_total
            qual_assess_df['FEATURE_ALIGN_NOT_UNIQUE'] = qual_assess_df['FEATURE_ALIGN_NOT_UNIQUE'] / unique_alignment_total
            qual_assess_df['AMBIGUOUS_FEATURE'] = qual_assess_df['AMBIGUOUS_FEATURE'] / unique_alignment_total
            qual_assess_df['TOO_LOW_AQUAL'] = qual_assess_df['TOO_LOW_AQUAL'] / unique_alignment_total

        try:
            if self.coverage_check_flag:
                coverage_df = self.coverageCheck()
                qual_assess_df = pd.merge(qual_assess_df, coverage_df, how='left', on='FASTQFILENAME')
        except AttributeError:
            self.logger.info('query_df or coverage_check_flag not present -- no coverage check')

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
            htseq_count_metadata_category = line[0][2:].upper()  # drop the __ in front of the category
            # enter to htseq_count_dict
            library_metadata_dict.setdefault(htseq_count_metadata_category, int(line[1]))
            # iterate
            line = next(htseq_file_reversed)
        # rename some key/value pairs
        library_metadata_dict['NOT_ALIGNED_TOTAL'] = library_metadata_dict.pop('NOT_ALIGNED')
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

    def coverageCheck(self):
        """
           calculate gene coverage of genes in the 'genotype' column of the query_df that do not have suffix _over
           split genotype into two columns, genotype1, genotype2 to address double KO
        """
        genotype_df = self.query_df[['fastqFileName', 'genotype']]
        genotype_df['fastqFileName'] = genotype_df['fastqFileName'].apply(lambda x: utils.pathBaseName(x))
        # new column perturbation, if _over in genotype, put 'over' otherwise 'ko'
        genotype_df['pertubation'] = ['over' if '_over' in genotype_column
                                        else 'ko' for genotype_column in genotype_df['genotype']]
        # remove _over from genotype
        genotype_df['genotype'] = genotype_df['genotype'].str.replace('_over', '')
        # split genotype on period. rename column 2 genotype2 if exists. if not, add genotype_2 with values None
        genotype_columns = genotype_df['genotype'].str.split('.', expand=True)
        if len(list(genotype_columns.columns)) == 2:
            genotype_columns.rename(columns={0:'genotype_1', 1: 'genotype_2'},inplace=True)
        else:
            genotype_columns.rename(columns={0:'genotype_1'},inplace=True)
            genotype_columns['genotype_2'] = None
        # bind genotype_df to genotype_columns
        genotype_df.drop(columns=['genotype'],inplace=True)
        genotype_df = pd.concat([genotype_df, genotype_columns], axis=1)
        # get bam filepaths
        try:  # TODO: clean this up. ugly method of testing whether a list of filepaths has been passed or not. this is for the nextflow pipeline
            bam_file_paths = [bam_file for bam_file in self.nextflow_list_of_files if '_sorted_aligned_reads.bam' in bam_file]
        except AttributeError:
            # extract files in directory with given suffix
            bam_file_paths = glob("{}/*{}".format(self.quality_assess_dir_path, '_sorted_aligned_reads.bam'))
        # set genomes
        s288c_r64_genome = os.path.join(self.genome_files, 'S288C_R64', 'S288C_R64.gff')
        if not os.path.isfile(s288c_r64_genome):
            self.logger.critical('S288C_R64 genome path not valid: %s' %s288c_r64_genome)
            raise FileNotFoundError
        kn99_genome = os.path.join(self.genome_files, 'KN99', 'crNeoKN99.gtf')
        if not os.path.isfile(kn99_genome):
            self.logger.critical('S288C_R64 genome path not valid: %s' %kn99_genome)
            raise FileNotFoundError
        h99_genome = os.path.join(self.genome_files, 'H99', 'crNeoH99.gtf')
        if not os.path.isfile(h99_genome):
            self.logger.critical('S288C_R64 genome path not valid: %s' %h99_genome)
            raise FileNotFoundError
        # create columns genotype_1_coverage and genotype_2_coverage
        genotype_df['genotype_1_coverage'] = None
        genotype_df['genotype_2_coverage'] = None
        # loop over rows, calculating coverage for each genotype (testing wither genotype2 is none and perturbation is _over
        for index, row in genotype_df.iterrows():
            if not row['genotype_1'] == 'CNAG_00000':
                # simple name is like this: run_673_s_4_withindex_sequence_TGAGGTT (no containing directories, no extention)
                print('...coverage check %s, %s' %(row['genotype_1'], row['genotype_2']))
                fastq_simple_name = utils.pathBaseName(row['fastqFileName'])
                genotype_1 = row['genotype_1']
                genotype_2 = row['genotype_2']
                if genotype_1.startswith('CNAG'):
                    genome = kn99_genome
                else:
                    genome = s288c_r64_genome
                try:
                    bam_file = [bam_file for bam_file in bam_file_paths if fastq_simple_name in bam_file][0]
                except IndexError:
                    self.logger.critical('bam file not found for %s' %fastq_simple_name)
                # cmd = "grep %s %s | grep CDS | gff2bed | samtools depth -a -b - %s | wc -l" %(genotype_1, genome, bam_file)
                # num_bases_in_cds = int(subprocess.getoutput(cmd))
                cmd = "grep %s %s | grep CDS | gff2bed | samtools depth -a -b - %s | grep -v 0 | wc -l" % (genotype_1, genome, bam_file)
                # num_bases_in_cds_with_one_or_more_read = int(subprocess.getoutput(cmd))
                # genotype_df.loc[index, 'genotype_1_coverage'] = num_bases_in_cds / float(num_bases_in_cds_with_one_or_more_read)
                print(cmd)
                if genotype_2 is not None:
                    cmd = "grep %s %s | grep CDS | gff2bed | samtools depth -a -b - %s | wc -l" % (genotype_2, genome, bam_file)
                    # num_bases_in_cds = int(subprocess.getoutput(cmd))
                    # cmd = "grep %s %s | grep CDS | gff2bed | samtools depth -a -b - %s | grep -v 0 | wc -l" % (genotype_2, genome, bam_file)
                    # num_bases_in_cds_with_one_or_more_read = int(subprocess.getoutput(cmd))
                    # genotype_df.loc[index, 'genotype_2_coverage'] = num_bases_in_cds / float(num_bases_in_cds_with_one_or_more_read)
            # set as attribute self.coverage_check
        genotype_df.columns = [column_name.upper() for column_name in genotype_df.columns]
        return genotype_df[['FASTQFILENAME', 'GENOTYPE_1_COVERAGE', 'GENOTYPE_2_COVERAGE']]

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
