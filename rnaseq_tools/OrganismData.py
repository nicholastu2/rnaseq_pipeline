from rnaseq_tools import utils
from rnaseq_tools import StandardData

class OrganismData(StandardData.StandardData):
    def __init__(self, **kwargs):
        # additional attributes to add to the _attributes in StandardData
        self._organism_attributes = ['organism_genome_files']
        # initialize Standard data with the extended _attributes
        super(OrganismData, self).__init__(self._organism_attributes, **kwargs)

