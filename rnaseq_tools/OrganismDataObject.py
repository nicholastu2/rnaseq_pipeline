from rnaseq_tools import utils
from rnaseq_tools.StandardDataObject import StandardData
import os


class OrganismData(StandardData):
    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = ['organism', 'output_dir', 'wildtype', 'experiment_dir', 'norm_count_path',
                                         'max_replicates', 'drug_marker', 'qc_config', 'experiment_conditions']
        # TODO: This is a messy and repetitive way of adding expected attributes from children of OrganismData to add to StandardData
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # set list of known organisms
        self._configured_organisms_list = ['H99', 'KN99', 'S288C_R64']
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(OrganismData, self).__init__(self._add_expected_attributes, **kwargs)

        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'OrganismData'

        # set organism, if an organism is passed
        if hasattr(self, 'organism'):
            if self.organism in self._configured_organisms_list:
                self.setOrganismData()
            else:
                print('\n{self.organism} is not configured. You will have to set the OrganismData attributes manually. '
                      'See the config/rnaseq_pipeline_config.ini. Alternatively, see one of the configured genome_files (in {self.genome_files}) '
                      'and create a subdir of genomes_files with an OrganismData_config.ini file, zip it into '
                      '/lts/mblab/Crypto/rnaseq_data/genome_files.zip, remove your genome_files in your {self.user_rnaseq_pipeline_directory}'
                      'and either re-run this script or start an interactive python session, import and instantiate a StandardData object.\n')
                # see [OrganismData] in config/rnaseq_pipeline_config.ini
                utils.configure(self, self.config_file, self.self_type)

        # create OrganismData logger
        try:
            self.logger = None
            self.createOrganismDataLogger()
        except NotADirectoryError:
            print('Cannot create OrganismData logger. Check code and config file.')
            exit(1)
        else:
            if not os.path.isfile(self.logger_file_path):
                raise FileNotFoundError('LoggerNotSuccessfullyCreated')

    # noinspection PyAttributeOutsideInit
    def setOrganismData(self):
        # first, run standard directory structure to check that file structure exists, attributes set, etc
        self.standardDirectoryStructure()
        # set OrganismData config found in rnaseq_pipeline/genome_files/<organism>/OrganismData_config.ini
        setattr(self, 'organism_config_file', os.path.join(self.user_rnaseq_pipeline_directory, self.genome_files,
                                                           self.organism, 'OrganismData_config.ini'))
        # read configuration file at path stored in organism_config_file
        utils.configure(self, self.organism_config_file, self.self_type, os.path.join(self.genome_files,
                                                                                      self.organism))
        # remove '' after taking basename (see todo regarding basename)
        self.feature_type = utils.pathBaseName(self.feature_type).replace('\'', '')  # TODO: fix!! this is a problem. All other features are being set to paths, but not this one. this is an issue with using utils.configure, it seems

    def createOrganismDataLogger(self):
        """
            create logger for OrganismData
            :raises: NotADirectoryError if logger_directory_path does not exist
        """
        logger_directory_path = utils.dirPath(self.logger_file_path)
        if os.path.isdir(logger_directory_path):
            utils.createLogger(self.logger_file_path, __name__)
        else:
            raise NotADirectoryError('LogDirectoryDoesNotExist')
