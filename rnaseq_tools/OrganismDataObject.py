from rnaseq_tools import utils
from rnaseq_tools.StandardDataObject import StandardData
import pandas as pd
import numpy as np
import os
import configparser


class OrganismData(StandardData):
    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = ['organism', 'output_dir', 'wildtype', 'experiment_dir', 'norm_count_path',
                                         'max_replicates', 'drug_marker', 'qc_config', 'experiment_conditions']
        # TODO: This is a messy and repetitive way of adding expected attributes from children of OrganismData to add to StandardData
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(OrganismData, self).__init__(self._add_expected_attributes, **kwargs)
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'OrganismData'

        # set organism, if an organism is passed
        if hasattr(self, 'organism'):
            # set organism directory
            self.organism_directory = os.path.join(self.user_rnaseq_pipeline_directory, self.genome_files,
                                                   self.organism)
            # set OrganismData config found in rnaseq_pipeline/genome_files/<organism>/OrganismData_config.ini
            self.organism_config_file = os.path.join(self.organism_directory, 'OrganismData_config.ini')
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
        self.logger = utils.createStandardObjectChildLogger(self, __name__)

    def setOrganismData(self):
        # read configuration file at path stored in organism_config_file
        # read config file
        config = configparser.ConfigParser()
        config.read(self.organism_config_file)
        # set attributes for StandardData
        for key, value in config[self.self_type].items():
            if key in self._no_file_organism_attributes['strings']:
                setattr(self, key, str(value))
            elif key in self._no_file_organism_attributes['ints']:
                setattr(self, key, int(value))
            else:
                setattr(self, key, os.path.join(self.organism_directory, value))

    def createOrganismDataLogger(self):
        """
            create logger for OrganismData
            :raises: NotADirectoryError if logger_directory_path does not exist
        """
        logger_directory_path = utils.dirPath(self.log_file_path)
        if os.path.isdir(logger_directory_path):
            self.logger = utils.createStandardObjectChildLogger(self, __name__)
        else:
            raise NotADirectoryError('LogDirectoryDoesNotExist')

    def createCountSheet(self, count_file_list):
        """
            create count matrix with list of genes as rows and samples as columns
            :param exp_dir: experiment directory created by create_experiment.py
            :param gene_list: see the OrganismData_config.ini in genome_files/<organism>
            :returns: A count matrix of all genes (rows) by all samples (columns)
        """
        # instantiate count_df with columns gene_id, <sample_name_1...>, <sample_name_2...>
        count_df = pd.DataFrame(columns=['gene_id'])
        # add gene names to column gene_id
        with open(self.gene_list) as gene_file:
            gene_list_generator = gene_file.readlines()
            gene_list = [gene_name.rstrip() for gene_name in gene_list_generator]
        # add gene names to gene_id
        count_df.gene_id = gene_list

        # ensure that items in column_list are only basenames, not paths
        column_list = [os.path.basename(x) for x in count_file_list]
        for i in range(len(count_file_list)):  # note: column_list is list of basepaths, count_file_list is list of full paths. need count_file_list item to read in file, but column names are in column_list
            # append sample column counts
            print('...working on %s' % column_list[i])
            count_file_df = pd.read_csv(count_file_list[i], sep='\t', names=['gene_id', column_list[i]])
            count_df = pd.merge(count_df, count_file_df, on='gene_id')

        return count_df
