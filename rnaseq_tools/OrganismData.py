from rnaseq_tools import utils
from rnaseq_tools.StandardData import StandardData
import configparser
import os

class OrganismData(StandardData):
    def __init__(self, **kwargs):
        # set (overwrite?) self_type depending on cmdline input
        self.self_type = 'OrganismData'
        # and configure
        utils.configure(self)
        self.list_of_known_organisms = ['H99', 'KN99', 'S288C_R64']
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(OrganismData, self).__init__(**kwargs)
        # make sure self_type wasn't overwritten by super call above
        self.self_type = 'OrganismData'
        # set organism data, if it is passed
        if hasattr(self, 'organism'):
            if self.organism in self.list_of_known_organisms:
                self.setOrganismData()
            else:
                print('\n{self.organism} is not configured. You will have to set the OrganismData attributes manually. '
                      'See the config/rnaseq_pipeline_config.ini. Alternatively, see one of the configured genome_files (in {self.genome_files}) '
                      'and create a subdir of genomes_files with an OrganismData_config.ini file, zip it into '
                      '/lts/mblab/Crypto/rnaseq_data/genome_files.zip, remove your genome_files in your {self.user_rnaseq_pipeline} directory'
                      'and either re-run this script or start an interactive python session, import and instantiate a StandardData object.\n')

    def setOrganismData(self):
        setattr(self, 'organism_config_file', os.path.join(self.user_rnaseq_pipeline, self.genome_files,
                                                           self.organism, 'OrganismData_config.ini'))
        utils.configure(self, self.organism_config_file, self.self_type)




