from rnaseq_tools import utils
from rnaseq_tools import StandardData

class IgvObject(StandardData.StandardData):
    def __init__(self, **kwargs):
        # additional attributes to add to the _attributes in StandardData
        self._igv_attributes = ['sample_list', 'igv_genome', 'output_dir']
        # initialize Standard data with the extended _attributes
        super(IgvObject, self).__init__(self._igv_attributes, **kwargs)



    def makeIgvDict(self):
        """
        from list of samples, create dictionary in format {'sample': {'gene': [gene_1, gene_2,...], 'bed': /path/to/bed, 'bam': /path/to/bam}
        """
        if not hasattr(self, self.sample_list):
            raise AttributeError('Your instance of IgvObject does not have a list of samlpes. '
                                 'It must have a list of samples (list of fastq or count file names) to create the igv_dict.')






def writeIgvJobScript(ineffmut_dict, igv_genome, igv_output_dir, fig_format='png', email=None,
                     job_script='job_scripts/igv_snapshot.sbatch'):
    """
	Write sbatch job script to make IGV snapshot
	"""
    num_samples = len(ineffmut_dict.keys())
    job = '#!/bin/bash\n#SBATCH -N 1\n#SBATCH --mem=5G\n'
    job += '#SBATCH -D ./\n#SBATCH -o log/igv_snapshot_%A.out\n#SBATCH ' \
           '-e log/igv_snapshot_%A.err\n#SBATCH -J igv_snapshot\n'
    if email is not None:
        job += '#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=%s\n' % email
    job += '\nml java\n'
    for sample in ineffmut_dict.keys():
        bam_file = ineffmut_dict[sample]['bam']
        bed_file = ineffmut_dict[sample]['bed']
        # this is a call to another script in the rnaseq_pipeline/tools
        job += 'python -u tools/make_IGV_snapshots.py %s ' \
               '-bin /opt/apps/igv/2.4.7/igv.jar -nf4 -r %s -g %s -fig_format %s -o %s\n' \
               % (bam_file, bed_file, igv_genome, fig_format, igv_output_dir)
    # write job to script
    writer = open(job_script, 'w')
    writer.write('%s' % job)
    writer.close()