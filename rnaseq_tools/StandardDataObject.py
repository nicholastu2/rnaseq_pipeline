from rnaseq_tools import utils
import os
import getpass  # see https://www.saltycrane.com/blog/2011/11/how-get-username-home-directory-and-hostname-python/


class StandardData:
    """
    parent class of rnaseq_pipeline. Creates rnaseq_pipeline directory in $USER if it (or any part) do not exist
    and stores all paths to resource directories and files.
       Does not check contents of all required directories for completeness (especially critical for genome_files), #TODO fix this so it does
       only that they exist.
       loads rnaseq_pipeline package level configuration (both main config as well as logger config in $USER/rnaseq_pipeline/config
           Recall that OrganismData, a child of StandardData, loads the .ini config files in each organism in $USER/rnaseq_pipeline/genome_files/
    """

    def __init__(self, expected_attributes=None, *args, **kwargs):
        """
        initialize StandardDataFormat with arbitrary number of keyword arguments.
        :param expected_attributes: a list of other attributes (intended use is for sub classes to add expected attributes)
        :param kwargs: arbitrary number/length keyword arguments. key = value will be set as class attributes
        """
        self.self_type = 'StandardData'
        self.log_dir = None
        self.log_file_path = None
        # list of StandardData object expected attributes. This is the file structure necessary for the rnaseq pipeline
        self._attributes = ['lts_rnaseq_data', 'pipeline_version', 'mblab_scratch', 'scratch_database_files',
                            'mblab_shared', 'lts_sequence', 'lts_align_expr', 'scratch_sequence',
                            'user_rnaseq_pipeline_directory',
                            'genome_files', 'reports', 'sbatch_log', 'log_dir', 'log_file', 'job_scripts', 'rnaseq_tmp',
                            'config_file', 'align_count_path']

        # these run numbers have leading zeros in either/or the database or lts_align_expr. See -lz in create_experiment.py input options
        self._run_numbers_with_zeros = {641: '0641', 647: '0647', 648: '0648', 659: '0659', 673: '0673', 674: '0674',
                                        684: '0684', 731: '0731', 748: '0748', 759: '0759', 769: '0769', 773: '0773',
                                        779: '0779'}

        # set year_month_day
        self.year_month_day = utils.yearMonthDay()

        # This is to extend _attributes if a class extends StandardData
        if isinstance(expected_attributes, list):
            self._attributes.extend(expected_attributes)

        # get user name and set as _user
        self._user = getpass.getuser()

        # set config, either entered when creating StandardDataObject or child, or use default expected on cluster # TODO: possible? to containerize scripts calling this with a repo structure and have path always point there (align_counts eg)
        self.default_config_path = '/opt/apps/labs/mblab/software/rnaseq_pipeline/1.0/config/rnaseq_pipeline_config.ini'
        try:
            self.config_file = kwargs['config_file']
            if not os.path.isfile(self.config_file):
                raise FileNotFoundError('ConfigFileNotFound')
        except KeyError or FileNotFoundError:
            try:
                if not os.path.exists(self.default_config_path):
                    raise FileNotFoundError('DefaultConfigPathNotValid')
            except FileNotFoundError:
                print('Either specify, or check the path to, config_file = /path/to/config/file in your call to StandardDataObject or Child')
        self.config_file = self.default_config_path
        utils.setAttributes(self, self._attributes, kwargs)
        # load config file
        utils.configure(self, self.config_file, self.self_type)

        # set interactive (flag for interactive session on htcf) to false if not already set. If True, StandardDataObject and child will not try to softlink to lts (long term storage)
        if not hasattr(self, 'interactive'):
            self.interactive = False  # TODO: NEED TO CHECK THIS -- CHILDREN MAY NOT OVERWITE!

    def standardDirectoryStructure(self):
        """
        checks for and creates if necessary the expected directory structure in /scratch/mblab/$USER/rnaseq_pipeline
        """
        # set attribute user_scratch (this is where rnaseq_pipeline and all subordinate folders/files will be
        user_scratch = os.path.join(self.mblab_scratch, self._user)
        try:
            if not os.path.isdir(user_scratch):
                raise NotADirectoryError
        except NotADirectoryError:
            print('%s is not a valid path to a directory.')
        else:
            setattr(self, 'user_scratch', user_scratch)

        # if it does not already exist, create user_rnaseq_pipeline in user_scratch and set attribute
        utils.mkdirp(self.user_rnaseq_pipeline_directory)
        setattr(self, 'user_rnaseq_pipeline_directory', '{}/rnaseq_pipeline'.format(self.user_scratch))

        # create necessary subdirectories in rnaseq_pipeline
        process_directories = ['reports', 'query', 'sbatch_log', 'log/%s' % self.year_month_day, 'job_scripts',
                               'rnaseq_tmp', 'experiments']  # TODO: MAKE SBATCH_LOG LIKE LOG WITH YEAR_MONTH_DAY SUBDIR
        for directory in process_directories:
            # store path
            path = os.path.join(self.user_rnaseq_pipeline_directory, directory)
            # this will only create the path if it dne
            utils.mkdirp(path)
            # set attr to directory (the names in process_directories) unless log, which is treated specially
            if directory == 'log/%s' % self.year_month_day:
                # distinguish the log directory ($USER/rnaseq_pipeline/log)
                self.log_dir = 'log/%s' % self.year_month_day
                utils.mkdirp(self.log_dir)
                # from the daily log file ($USER/rnaseq_pipeline/log/<year-month-day>)
                self.log_file_path = os.path.join(self.log_dir, '%s.log' % self.year_month_day)
            else:
                setattr(self, directory, path)

        # TODO: priority (see other comment)
        # check for the directories to be soft linked from /scratch/mblab/mblab.shared (self.mblab_shared). soft link if not, setattr in either case # TODO: decouple this from scratch space -- maybe pull from database-files repo and from path to genome_files rather than creating that from lts
        try:
            mblab_shared_dirs = ['scratch_sequence', 'database_files']
            utils.softLinkAndSetAttr(self, mblab_shared_dirs, self.mblab_shared, self.user_rnaseq_pipeline_directory)
        except FileNotFoundError:
            print("Could not softlink scratch_sequence or database_files. check both source and target.\n"
                  "To check target, look in StandardDataObject.standardDirectoryStructure() and config_file.")

        # if interactive is False, link from lts # TODO: figure out how to do this while off the cluster -- change interactive to align_expr_sequence or something like that, set false if don't have access to those directories
        if not self.interactive:
            # check for directories to be soft linked from /lts/mblab/Crypto/rnaseq_pipeline (self.lts_rnaseq_data)
            lts_dirs_to_softlink = ['lts_align_expr', 'lts_sequence']
            utils.softLinkAndSetAttr(self, lts_dirs_to_softlink, self.lts_rnaseq_data,
                                     self.user_rnaseq_pipeline_directory)

        # TODO: priority figure out how to do this without pulling from /lts. put link to genome_files.zip in config maybe
        # unzip genome files from /lts/mblab/Crypto/rnaseq_data/1.0/genome_files to self.user_rnaseq_pipeline_directory
        setattr(self, 'genome_files', os.path.join(self.user_rnaseq_pipeline_directory, 'genome_files'))
        if not os.path.exists(self.genome_files):
            genome_files_full_path = os.path.join(self.lts_rnaseq_data, self.pipeline_version, 'genome_files.zip')
            cmd = 'unzip {} -d {}'.format(genome_files_full_path, self.user_rnaseq_pipeline_directory)
            utils.executeSubProcess(cmd)

    def setStandardDataLogger(self):
        """
            function to create StandardData logger
        """
        try:
            setattr(self, 'logger', utils.createLogger(self.log_file_path, __name__))
        except NameError:
            print('cannot set logger without specifying log_file_path in StandardDataObject/child and self.standardDirectoryStructure()')
            exit(1)
