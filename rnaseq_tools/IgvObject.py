"""
tools to create browser shots

makeIgvSnapshotDict create a dict with the following format
igv_snapshot_dict = {'wildtype.fastq.gz': {'gene': [gene_1, gene_2, gene_3...],
                     'bam': file.bam, 'bed': file.bed},
                     'perturbed_sample.fastq.gz': {'gene': [gene_1, gene_2],
                     'bam': file.bam, 'bed': file.bed},
                     'perturbed_sample.fastq.gz': {'gene': [gene_1],
                     'bam': file.bam, 'bed': file.bed}
                  }
and ensures that all alignment files (.bam) and their index companions (.bam.bai) are in the experiment directory, which
is either supplied or created in rnaseq_tmp

usage: igv = IgvObject(organism = <your organism>, query_sheet_path = <path to query>,
                       sample_list = a python list of samples to take shot of,
                       igv_output_dir = output directory for igv snapshots,
                       experiment_dir = <path to experiment dir> # if this is not present, bed/bat files will be put in rnaseq_tmp/datetime_igv_files)

"""

from rnaseq_tools import utils
from rnaseq_tools.OrganismData import OrganismData
from rnaseq_tools import annotation_tools
import sys
import os
import time
import subprocess


# igv_snapshot_dict = {'wildtype.fastq.gz': {'gene': [gene_1, gene_2, gene_3...], 'bam':
#                  file.bam, 'bed': file_1.bed},
#                  'perturbed_sample.fastq.gz': {'gene': [gene_1, gene_2], 'bam':
#                   file.bam, 'bed': file_1.bed},
#                  'perturbed_sample.fastq.gz': {'gene': [gene_1], 'bam': file.bam,
#                  'bed': file_1.bed}
#                  }

class IgvObject(OrganismData):
    """
    :attribute _igv_attributes: list of expected attributes of a IgvObject. Like other StandardData children, these will update
    self._attribute (inherited from StandardData through any StandardData to current object). These attributes are attributes which
    some level of the StandardData class will recognize and know what to do with. They are variables that are expected to be either
    entered by a config file or by the user on instantiation. If they are not set, *typically* they will be set automatically in the script (this is a TODO OBJECT TO TAKE CARE OF THIS BETTER)
    or a print statement will tell you the default location, or an error will be raised (again, TODO)
    """

    def __init__(self, **kwargs):
        # additional attributes to add to the _attributes in StandardData
        # TODO: possibly change inheretence to a subclass of OrganismData that sets up a class for ANY scheduler manipulation (ie align_counts, this) that take email as an optional argument
        self._igv_attributes = ['sample_list', 'igv_genome', 'igv_output_dir']
        # initialize Standard data with the extended _attributes
        super(IgvObject, self).__init__(self._igv_attributes, **kwargs)
        # initialize list to store bamfiles that need to be indexed (must be done by batch script)
        self.bam_file_to_index_list = []

    def checkAttributes(self):
        if not hasattr(self, 'sample_list'):
            raise AttributeError('Your instance of IgvObject does not have a list of samples. '
                                 'It must have a list of samples (list of fastq or count file names) to create the igv_dict. '
                                 '\nThis may be just a list version of the COUNTFILENAME column of self.query_df')
        # raise attribute error if no (standardized) query_df. this should be created by passing query_sheet_path to the constructor. See StandardData
        if not hasattr(self, 'query_df'):
            raise AttributeError('No query_df (standardized query_df object in StandardData object). This is necessary '
                                 'and created if you input a query_sheet_path to the constructor of any StandardData '
                                 'object instance. However, you can input the full database (queryDB.py with flag -pf). '
                                 'It will just make the searches for a given sample somewhat longer (though likely not much).')
        # if an experiment_dir is not passed in constructor of IgvObject, create one in rnaseq_tmp/<timenow>_igv_files and store the path in self.experiment_dir
        if not hasattr(self, 'experiment_dir'):
            print(
                'creating a directory in rnaseq_tmp to store alignment files so that the scheduler has access to them.'
                'This assumes that the count and alignment files have already been processed and moved to '
                '/lts/mblab/Crypto/rnaseq_data/align_expr')
            timestr = time.strftime("%Y%m%d_%H%M%S")
            setattr(self, 'experiment_dir', os.path.join(self.rnaseq_tmp, '{}_igv_files'.format(timestr)))
            utils.mkdirp(self.experiment_dir)

    def moveAlignmentFiles(self): # TODO: THIS ALSO CHECKS FOR AND MAKES A LIST OF ALIGNMENT FILES TO INDEX WITH SAMTOOLS -- NOT CLEAR IN NAME
        """ TODO: a general move files script needs to be written. input list of files to be moved, source, dest, move files if they are not in dest. put in rnaseq_tools.utils
        ensure that all alignment files are in self.experiment_dir. If no self.experiment_dir set, then these will be deposited in
        rnaseq_tmp/datetime_igv_files
        This needs to be done b/c indexing can only take place via srun/sbatch (on compute node. samtools not available on login node as of 3/2020)
        """

        for sample in self.sample_list:
            # strip .fastq.gz if it is there and add _read_count.tsv -- remember that self.query_df which has fastqFileName --> COUNTFILENAME and .fastq.gz converted to _read_count.tsv extensions. All column headings CAPITAL
            if sample.endswith('.fastq.gz'):
                sample = utils.pathBaseName(sample) + '_read_count.tsv'
            # extract run_number just in case needed to find bam file in align_expr
            run_number = self.extractValueFromStandardizedQuery('COUNTFILENAME', sample, 'RUNNUMBER',
                                                                check_leading_zero=True)
            # create bamfile name
            bamfile = sample.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
            # if it is not in the exp dir, then add it
            if not os.path.exists(os.path.join(self.experiment_dir, bamfile)):
                prefix = utils.addForwardSlash(self.lts_align_expr)
                bamfile_align_expr_path = '{}run_{}/{}'.format(prefix, run_number, bamfile)
                cmd = 'rsync -aHv {} {}'.format(bamfile_align_expr_path, utils.addForwardSlash(self.experiment_dir))
                utils.executeSubProcess(cmd)
            # store full path to bam file in experiment dir
            bamfile_fullpath = os.path.join(self.experiment_dir, bamfile)
            if not os.path.exists(bamfile_fullpath + '.bai'):  # test if indexed bam exists
                self.bam_file_to_index_list.append(bamfile_fullpath)

    def writeIndexScript(self):
        """
       write sbatch script to index bam files
       This needs to be done b/c indexing can only take place via srun/sbatch (on compute node. samtools not available on login node as of 3/2020)
       script will be place in rnaseq_pipeline/job_scripts

       """
        self.igv_index_script = os.path.join(self.job_scripts, 'index_alignment_files.sbatch')

        job = '#!/bin/bash\n' \
              '#SBATCH -N 1\n' \
              '#SBATCH --mem=5G\n' \
              '#SBATCH -o {0}/index_bams_%A.out\n' \
              '#SBATCH -e {0}/index_bams_%A.err\n' \
              '#SBATCH -J index_bams\n'.format(self.log)
        if hasattr(self, 'email'):
            job += '#SBATCH --mail-type=END,FAIL\n' \
                   '#SBATCH --mail-user=%s\n' % self.email
        job += '\nml samtools\n'

        for alignment_file in self.bam_file_to_index_list:
            job += '\nsamtools index -b %s -o %s\n' %(alignment_file, self.experiment_dir)

        with open(self.igv_index_script, 'w') as file:
            file.write('%s' % job)


    def makeIgvSnapshotDict(self):
        """
        from list of samples, create dictionary in format {'sample': {'gene': [gene_1, gene_2,...], 'bed': /path/to/bed, 'bam': /path/to/bam}
        # ASSUMPTION: THE SAMPLE LIST IS EITHER A LIST OF .FASTQ.GZ OR *_read_count.tsv
        """
        igv_snapshot_dict = {}
        genotype_list = []
        wildtype_sample_list = []
        for sample in self.sample_list:
            bamfile = sample.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
            bamfile_fullpath = os.path.join(self.experiment_dir, bamfile)
            if not os.path.exists(bamfile_fullpath):
                print('bamfile does not exist in experiment_directory %s. Run moveAlignmentFiles first.')
                break
            genotype = self.extractValueFromStandardizedQuery('COUNTFILENAME', sample, 'GENOTYPE')
            # split on period if this is a double perturbation. Regardless of whether a '.' is present,
            # genotype will be cast to a list eg ['CNAG_00000'] or ['CNAG_05420', 'CNAG_01438']
            genotype = genotype.split('.')
            # if genotype ends with _over, remove _over
            for index in range(len(genotype)):
                genotype[index] = genotype[index].replace('_over', '')
            # if the object has an attribute wildtype, and genotype is not wildtype, add to igv_snapshot_dict
            if hasattr(self, 'wildtype') and genotype[0] == self.wildtype:
                # add the genotypes to genotype_list (not as a list of list, but as a list of genotypes)
                genotype_list.extend(genotype)
                # add to igv_snapshot_dict
                igv_snapshot_dict.setdefault(sample, {}).setdefault('gene', []).extend(
                    genotype)  # TODO: clean this up into a single line, make function to avoid repeated code below w/wildtype
                igv_snapshot_dict[sample]['bam'] = bamfile_fullpath
                igv_snapshot_dict[sample]['bed'] = None
            # if genotype is equal to wildtype, then store the sample as the wildtype (only one, check if this is right)
            else:
                igv_snapshot_dict.setdefault(sample, {'gene': None, 'bam': None, 'bed': None})
                wildtype_sample_list.append([sample, bamfile_fullpath])  # wildtype_sample_list will be a list of lists
        # if the wildtype genotype was found, create entry in the following form
        # {sample_read_counts.tsv: {'gene': [perturbed_gene_1, perturbed_gene_2, ...], 'bam': wt.bam, 'bed': created_bed.bed}
        for wt_sample in wildtype_sample_list:
            igv_snapshot_dict[wt_sample[0]]['gene'] = genotype_list
            igv_snapshot_dict[wt_sample[0]]['bam'] = wt_sample[1]
            igv_snapshot_dict[wt_sample[0]]['bed'] = None
        # set attribute pointing toward igv_snapshot_dict
        setattr(self, 'igv_snapshot_dict', igv_snapshot_dict)

    # def indexBam(self, bamfile_fullpath):
    #     """
    #     Index the bam files. The .bai file (indexed bam) will be deposited in the same directory as the bam file itself.
    #     for igv, as long as the .bai is in the same directory as the .bam file, it will work.
    #     :param bamfile_fullpath: full path (absolute) to bamfile
    #     """
    #     print('\nindexing {}'.format(bamfile_fullpath))
    #     # subprocess.call will wait until subprocess is complete
    #     exit_status = subprocess.call('samtools index -b {}'.format(bamfile_fullpath), shell=True)
    #     if exit_status == 0:
    #         print('\nindexing complete. .bai is deposited in {}'.format(self.experiment_dir))
    #     else:
    #         sys.exit('failed to index {}. Cannot continue'.format(bamfile_fullpath))

    def createBedFile(self, flanking_region=500, file_format='png'):
        """
        Create bed files to describe IGV region of interest
        :param flanking_region: how far up/down stream from the gene to capture in the snapshot
        :param file_format: what format to use for image file
        """
        # TODO if bed files exist in exp dir, just get those

        ## get gene dictionary with chromsome, gene coordinates, strand
        if self.annotation_file.endswith('gtf'):
            self.annotation_dict = annotation_tools.parseGtf(self.annotation_file)
        elif self.annotation_file.endswith('gff') or self.annotation_file.endswith('gff3'):
            self.annotation_dict = annotation_tools.parseGff3(self.annotation_file)
        else:
            sys.exit(
                "ERROR: The gene annotation format cannot be recognized.")  # TODO: clean up preceeding blocks -- move parseGFF to OrganismData
        ## create gene body region bed file
        for sample in self.igv_snapshot_dict.keys():
            igv_bed_filepath = os.path.join(self.experiment_dir, utils.pathBaseName(sample) + '.bed')
            self.igv_snapshot_dict[sample]['bed'] = igv_bed_filepath
            with open(igv_bed_filepath, 'w') as file:
                for gene in self.igv_snapshot_dict[sample]['gene']:
                    d = self.annotation_dict[gene]
                    file.write('%s\t%d\t%d\t[%s]%s.%s\n' % (d['chrm'], d['coords'][0] - flanking_region,
                                                            d['coords'][1] + flanking_region, sample, gene,
                                                            file_format))

    def writeIgvJobScript(self, fig_format='png'):
        """
        Write sbatch job script to make IGV snapshot
        """
        # set igv_output_dir if not entered on cmd line
        if not hasattr(self, 'igv_output_dir'):
            setattr(self, 'igv_output_dir', os.path.join(self.experiment_dir, 'igv_output'))
        num_samples = len(self.igv_snapshot_dict.keys())
        job = '#!/bin/bash\n' \
              '#SBATCH -N 1\n' \
              '#SBATCH --mem=5G\n' \
              '#SBATCH -o {0}/igv_snapshot_%A.out\n' \
              '#SBATCH -e {0}/igv_snapshot_%A.err\n' \
              '#SBATCH -J igv_snapshot\n'.format(self.log)
        if hasattr(self, 'email'):
            job += '#SBATCH --mail-type=END,FAIL\n' \
                   '#SBATCH --mail-user=%s\n' % self.email
        job += '\nml java\n' \
               'ml rnaseq_pipeline\n'
        for sample in self.igv_snapshot_dict.keys():
            bam_file = self.igv_snapshot_dict[sample]['bam']
            bed_file = self.igv_snapshot_dict[sample]['bed']
            # this is a call to another script in the rnaseq_pipeline/tools
            job += '\nmake_IGV_snapshots.py %s -bin /opt/apps/igv/2.4.7/igv.jar -nf4 -r %s -g %s -fig_format %s -o %s\n' \
                   % (bam_file, bed_file, self.genome_files, fig_format, self.igv_output_dir)
        # write job to script
        igv_job_script_path = os.path.join(self.job_scripts, utils.pathBaseName(self.experiment_dir) + '.sbatch')
        setattr(self, 'igv_job_script', igv_job_script_path)
        with open(self.igv_job_script, 'w') as file:
            file.write('%s' % job)
        # TODO add a log statement, think about redirecting these job scripts to the job_script directory
        print('the job scripts have been deposited in {}\n'
              'this has been set as the igv_output_directory'.format(self.igv_output_dir))

    # def takeSnapshot(self):
    #     # TODO: check that attr exist
    #
    #     # TODO: with subprocess.call, can we pass python data? If so, don't write this to file
    #     batchscript_file = os.path.join(outdir, "IGV_snapshots.bat")
    #     # TODO: this should be to log, not to stdout
    #     print('\n~~~ IGV SNAPSHOT AUTOMATOR ~~~\n')
    #     print('Reference genome:\n{}\n'.format(genome))
    #     print('Track height:\n{}\n'.format(image_height))
    #     print('IGV binary file:\n{}\n'.format(igv_jar_bin))
    #     print('Output directory will be:\n{}\n'.format(outdir))
    #     print('Batchscript file will be:\n{}\n'.format(batchscript_file))
    #     print('Region file:\n{}\n'.format(region_file))
    #     print('Input files to snapshot:\n')


x = IgvObject()
