from rnaseq_tools import utils
import os
import sys
import getpass  # see https://www.saltycrane.com/blog/2011/11/how-get-username-home-directory-and-hostname-python/
import configparser

class StandardData:
    """
        parent class of rnaseq_pipeline. Creates rnaseq_pipeline directory in $USER if it (or any part) do not exist
        and stores all paths to resource directories and files.
           loads rnaseq_pipeline package level configuration (both main config as well as logger config in $USER/rnaseq_pipeline/config
               Recall that OrganismData, a child of StandardData, loads the .ini config files in each organism in $USER/rnaseq_pipeline/genome_files/

        In the config file, if using locally, you can include genome_files = https://... path to genome files in /lts. See this:
            https://htcfdocs.readthedocs.io/en/latest/storage/#publishing-files
    """

    def __init__(self, expected_attributes=None, *args, **kwargs):
        """
            initialize StandardDataFormat with arbitrary number of keyword arguments.
            :param expected_attributes: a list of other attributes (intended use is for sub classes to add expected attributes)
            :param kwargs: arbitrary number/length keyword arguments. key = value will be set as class attributes
        """
        self.self_type = 'StandardData'
        # these will be filled in the standardDataStructure() call at the end of the constructor
        self.log_dir = None
        self.log_file_path = None

        # list of StandardData object expected attributes. This is the file structure necessary for the rnaseq pipeline
        self._attributes = ['lts_rnaseq_data', 'pipeline_version', 'mblab_scratch', 'scratch_database_files',
                            'mblab_shared', 'lts_sequence', 'lts_align_expr', 'scratch_sequence',
                            'user_rnaseq_pipeline_directory', 'genome_files', 'reports', 'align_count_results',
                            'sbatch_log', 'log_dir', 'log_file', 'job_scripts', 'rnaseq_tmp', 'config_file']

        # these run numbers have leading zeros in either/or the database or lts_align_expr. See -lz in create_experiment.py input options
        self._run_numbers_with_zeros = {641: '0641', 647: '0647', 648: '0648', 659: '0659', 673: '0673', 674: '0674',
                                        684: '0684', 731: '0731', 748: '0748', 759: '0759', 769: '0769', 773: '0773',
                                        779: '0779', 711: '0711_5_0718_7', 718: '0711_5_0718_7', 711507187: '0711_5_0718_7', 6290618: '0629_0618'}

        # list of organisms with configured subdirectory in genome_files
        self._configured_organisms_list = ['H99', 'KN99', 'S288C_R64']

        # used in checkGenomefiles() and in OrganismDataObject
        self._no_file_organism_attributes = {'strings': ['organism_genome_file', 'feature_type'],
                                             'ints': ['total_exon_bases', 'total_intergenic_bases', 'nat_cds_length', 'g418_cds_length']} # NOTE: reading in the config seems to cast these to lower?

        # set year_month_day
        self.year_month_day = utils.yearMonthDay()

        # This is to extend _attributes if a class extends StandardData
        if isinstance(expected_attributes, list):
            self._attributes.extend(expected_attributes)

        # get user name and set as _user
        self._user = getpass.getuser()

        # set debug level -- all SD objects may adjust by passing logger_level argument in constructor
        try:
            self.logger_level = kwargs['logger_level']
        except KeyError:
            self.logger_level = 'INFO'

        # set config, either entered when creating StandardDataObject or child, or use default expected on cluster # TODO: possible? to containerize scripts calling this with a repo structure and have path always point there (align_counts eg)
        self.default_config_path = '/opt/apps/labs/mblab/software/rnaseq_pipeline/1.0/config/rnaseq_pipeline_config.ini'
        # this is the default in scripts using StandardData. see qual_assess_1 for example
        try:
            if kwargs['config_file'] == '/see/standard/data/invalid/filepath/set/to/default':
                kwargs['config_file'] = self.default_config_path
        except KeyError:
            pass
        # set config file
        try:
            self.config_file = kwargs['config_file']
        except KeyError:
            self.config_file = self.default_config_path
        finally:
            if not os.path.isfile(self.config_file):
                sys.exit('Default path to the htcf config not valid. Either specify, or check the path to, config_file = /path/to/config/file in your call to StandardDataObject or Child')
            else:
                utils.setAttributes(self, kwargs)
                # load config file
                utils.configure(self, self.config_file, self.self_type)

        # set interactive (flag for interactive session on htcf) to false if not already set. If True, StandardDataObject and child will not try to softlink to lts (long term storage)
        if not hasattr(self, 'interactive'):
            self.interactive = False

        # set/check standardDirectoryStructure
        self.standardDirectoryStructure()

    def standardDirectoryStructure(self):
        """
            checks for and creates if necessary the expected directory structure in /scratch/mblab/$USER/rnaseq_pipeline
        """
        # offer method to set user_scratch in config file
        try:
            if not os.path.isdir(self.user_scratch):
                raise NotADirectoryError('UserScratchDirectoryNotPresent')
        except AttributeError:
            # set attribute user_scratch (this is where rnaseq_pipeline and all subordinate folders/files will be
            user_scratch = os.path.join(self.mblab_scratch, self._user)
            setattr(self, 'user_scratch', user_scratch)
        except NotADirectoryError:
            utils.mkdirp(self.user_scratch)

        # if it does not already exist, create user_rnaseq_pipeline in user_scratch and set attribute
        setattr(self, 'user_rnaseq_pipeline_directory', '{}/rnaseq_pipeline'.format(self.user_scratch))
        utils.mkdirp(self.user_rnaseq_pipeline_directory)

        # create necessary subdirectories in rnaseq_pipeline
        process_directories = ['reports', 'align_count_results', 'query', 'sbatch_log', 'log/%s' % self.year_month_day, 'job_scripts',
                               'rnaseq_tmp', 'experiments', 'scratch_sequence']  # TODO: MAKE SBATCH_LOG LIKE LOG WITH YEAR_MONTH_DAY SUBDIR
        for directory in process_directories:
            # store path
            path = os.path.join(self.user_rnaseq_pipeline_directory, directory)
            # this will only create the path if it dne
            utils.mkdirp(path)
            # set attr to directory (the names in process_directories) unless log, which is treated specially
            if directory == 'log/%s' % self.year_month_day:
                # distinguish the log directory ($USER/rnaseq_pipeline/log)
                self.log_dir = os.path.join(self.user_rnaseq_pipeline_directory, 'log/%s' % self.year_month_day)
                utils.mkdirp(self.log_dir)
                # from the daily log file ($USER/rnaseq_pipeline/log/<year-month-day>)
                self.log_file_path = os.path.join(self.log_dir, '%s.log' % self.year_month_day)
                self.createStandardDataLogger()
            else:
                setattr(self, directory, path)

        try:
            database_files_path = os.path.join(self.user_rnaseq_pipeline_directory, 'database_files')
            if not os.path.isdir(database_files_path):
                raise NotADirectoryError('DatabaseFilesNotFound: %s' %database_files_path)
        except NotADirectoryError:
            cmd = 'git clone https://github.com/BrentLab/database_files.git %s' %database_files_path
            utils.executeSubProcess(cmd)
        finally:
            setattr(self, 'database_files', database_files_path)

        if self.interactive:
            print('Remember you will not be able to access lts_align_expr or lts_sequence in an interactive session on htcf')
        else:
            # check for directories to be soft linked from /lts/mblab/Crypto/rnaseq_pipeline (self.lts_rnaseq_data)
            lts_dirs_to_softlink = ['lts_align_expr', 'lts_sequence']
            try:
                utils.softLinkAndSetAttr(self, lts_dirs_to_softlink, self.lts_rnaseq_data,
                                         self.user_rnaseq_pipeline_directory)
            except FileNotFoundError:
                print('WARNING: The source of %s does not exist and are not accessible. In the future, it is better to include the flag\n'
                      'interactive=True in the constructor of a StandardData object when you are in an interactive session.' %lts_dirs_to_softlink)
                setattr(self, 'lts_align_expr', os.path.join(self.user_rnaseq_pipeline_directory, 'lts_align_expr'))
                setattr(self, 'lts_sequence', os.path.join(self.user_rnaseq_pipeline_directory, 'lts_sequence'))
            # TODO: priority figure out how to do this without pulling from /lts. put link to genome_files.zip in config maybe

        # unzip genome files from /lts/mblab/Crypto/rnaseq_data/1.0/genome_files to self.user_rnaseq_pipeline_directory
        self.setGenomeFiles()
        # check that all files present in the OrganismDataConfig.ini file in the subdirectories of genome_files exist
        try:
            self.checkGenomeFiles()
        except NotADirectoryError:
            print('Genome Files are incomplete. Delete genome_files completely and re-run StandardDataObject or child '
                  'to re-download genome_files.\nNote: this cannot be done from an interactive session on HTCF.')
        except FileNotFoundError:
            print('Genome Files are incomplete. Delete genome_files completely and re-run StandardDataObject or child '
                  'to re-download genome_files.\nNote: this cannot be done from an interactive session on HTCF.')

    def setGenomeFiles(self):
        """
            set genome_files path and download, if genome_files  DNE. If config_file has genome_files = https://...
            Then the zip file will be downloaded from that path
            TODO: error checking if the config_file https path doesn't work
        """
        # if genome_files is set in config file
        if hasattr(self, 'genome_files'):
            # if the config_file has an entry genome_files = 'https://...' (link to the hosted genome files in /lts -- it is important that there be a single source for genome_files)
            if self.genome_files.startswith('https'):
                # and the file genome_files DNE in user_rnaseq_pipeline_directory, download from path
                if not os.path.isdir(os.path.join(self.user_rnaseq_pipeline_directory, 'genome_files')):
                    zipped_genome_files_path = os.path.join(self.user_rnaseq_pipeline_directory, 'genome_files.zip')
                    download_genome_files_cmd = 'wget -O %s %s' %(zipped_genome_files_path, self.genome_files)
                    utils.executeSubProcess(download_genome_files_cmd)
                    unzip_genome_files_cmd = 'unzip %s -d %s && rm %s' %(zipped_genome_files_path, self.user_rnaseq_pipeline_directory,
                                                                         zipped_genome_files_path)
                    utils.executeSubProcess(unzip_genome_files_cmd)

        # set path of self.genome_files to subdir of user_rnaseq_pipeline directory
        setattr(self, 'genome_files', os.path.join(self.user_rnaseq_pipeline_directory, 'genome_files'))

        # if the file DNE and interactive flag is set to False (not in interactive session on  htcf), then download from /lts
        if not (self.interactive or os.path.exists(self.genome_files)):
            genome_files_full_path = os.path.join(self.lts_rnaseq_data, self.pipeline_version, 'genome_files.zip')
            cmd = 'unzip {} -d {}'.format(genome_files_full_path, self.user_rnaseq_pipeline_directory)
            utils.executeSubProcess(cmd)

    def checkGenomeFiles(self):  # NOTE: need to update OrganismDataObject to expect this function TODO: IMPROVE ERROR CHECKING AND LOGGING. CLUNKY CURRENTLY.
        """
            read in OrganismDataConfig.ini from each expected subdir of genome_files and check if the path is valid.
            If it is not, ask user to check genome_files and/or delete genome_files and allow StandardDataObject
            to re-download to update paths
        """
        # list of attributes not to be checked for file existence as they are not files
        no_check_organism_attribute_list = self._no_file_organism_attributes.values()
        # flatten list
        no_check_organism_attribute_list = [x for sublist in no_check_organism_attribute_list for x in sublist]
        for organism in self._configured_organisms_list:
            # check if directory exists
            organism_genome_files_subdir_path = os.path.join(self.genome_files, organism)
            if not os.path.isdir(organism_genome_files_subdir_path):
                self.logger.warning('%s does not exist in genome_files' %organism_genome_files_subdir_path)
                raise NotADirectoryError('ConfiguredOrganismSubdirectoryNotPresentInGenomeFiles')
            # check if the organism config file exists
            organism_config_file_path = os.path.join(organism_genome_files_subdir_path, 'OrganismData_config.ini')
            if not os.path.isfile(organism_config_file_path):
                self.logger.warning('The OrganismData_config.ini file does not exist for %s' % organism)
                raise FileNotFoundError('OrganismDataConfigFileNotFound')
            # if it does, read it in
            else:
                # read in config file as dictionary {genome_file_attribute: filename, ...} eg {novoalign_index: KN99_novoalign.nix}
                organism_config_dict = configparser.ConfigParser()
                organism_config_dict.read(organism_config_file_path)
                for organism_attribute, filename in organism_config_dict['OrganismData'].items():
                    # skip attributes that do not have a corresponding filename
                    if organism_attribute in no_check_organism_attribute_list:  # TODO: just make it so no int values are checked?
                        continue
                    organism_attribute_filepath = os.path.join(organism_genome_files_subdir_path, filename)
                    if not os.path.isfile(organism_attribute_filepath):
                        self.logger.warning('%s not found in %s subdirectory of genome_files' % (organism_attribute, organism))
                        raise FileNotFoundError('OrganismFileNotFound: %s for %s' %(organism_attribute, organism))

    def createStandardDataLogger(self):
        """
            function to create StandardData logger
        """
        try:
            setattr(self, 'logger', utils.createLogger(self.log_file_path, __name__, self.logger_level))
        except NameError:
            print('cannot set logger without specifying log_file_path in StandardDataObject/child and self.standardDirectoryStructure()')
            exit(1)

    def extractRunNumber(self, run_number):
        """
            input run number. if it is in the _run_numbers_with_leading_zero dict, then return a string with the leading
            zero
        """
        try:
            return self._run_numbers_with_zeros[int(run_number)] # how to handle casting better?
        except KeyError:
            return run_number
