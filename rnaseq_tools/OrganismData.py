from rnaseq_tools import utils
from rnaseq_tools import StandardData
import configparser

class OrganismData(StandardData.StandardData):
    def __init__(self, **kwargs):
        # set (overwrite?) self_type depending on cmdline input
        self.self_type = 'ObjectData'
        # and configure
        utils.configure(self)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(OrganismData, self).__init__(**kwargs)
        # make sure self_type wasn't overwritten by super call above
        self.self_type = 'ObjectData'
        # set organism data, if it is passed
        if hasattr(self, 'organism'):
            self.setOrganismData()

    def setOrganismData(self):
        config = configparser.ConfigParser()
        config.read(self.config_file)
        for organism_attribute in config['OrganismData']:
            # TODO grep the file from the file in genome files using a dictionary associating a organism_attribute with a file extension
            path = 'see above'
            setattr(self, organism_attribute, '/path/to/genome_files/{}/{}'.format(self.organism, path))




