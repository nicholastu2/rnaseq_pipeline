from rnaseq_tools import utils
from rnaseq_tools import StandardData
import pandas as pd
import pysam
import sys
import os

#in current code ineffmut_dict
# igv_snapshot_dict = {'wildtype.fastq.gz': {'gene': [gene_1, gene_2, gene_3...], 'bam':
#                  file.bam, 'bed': file_1.bed},
#                  'perturbed_sample.fastq.gz': {'gene': [gene_1, gene_2], 'bam':
#                   file.bam, 'bed': file_1.bed},
#                  'perturbed_sample.fastq.gz': {'gene': [gene_1], 'bam': file.bam,
#                  'bed': file_1.bed}
#                  }

class IgvObject(StandardData.StandardData):
    def __init__(self, **kwargs):
        # additional attributes to add to the _attributes in StandardData
        self._igv_attributes = ['sample_list', 'igv_genome', 'output_dir', 'wt', 'experiment_dir']
        # initialize Standard data with the extended _attributes
        super(IgvObject, self).__init__(self._igv_attributes, **kwargs)



    def makeIgvSnapshotDict(self):
        """
        from list of samples, create dictionary in format {'sample': {'gene': [gene_1, gene_2,...], 'bed': /path/to/bed, 'bam': /path/to/bam}
        # ASSUMPTION: THE SAMPLE LIST IS EITHER A LIST OF .FASTQ.GZ OR *_read_count.tsv
        """
        if not hasattr(self, self.sample_list):
            raise AttributeError('Your instance of IgvObject does not have a list of samlpes. '
                                 'It must have a list of samples (list of fastq or count file names) to create the igv_dict.')
        if not hasattr(self, self.query_df):
            raise AttributeError('No query_df (standardized query_df object in StandardData object). This is necessary '
                                 'and created if you input a query_sheet_path to the constructor of any StandardData '
                                 'object instance. However, you can input the full database (queryDB.py with flag -pf). '
                                 'It will just make the searches for a given sample somewhat longer (though likely not much).')
        igv_snapshot_dict = {}
        genotype_list = []
        bam_file_list = []
        for sample in self.sample_list:
            if sample.endswith('.fastq.gz'):
                sample = utils.pathBaseName(sample) + '_read_count.tsv'
            query_row_with_sample = super.query_df[super.query_df['COUNTFILENAME'] == sample]
            # extract value in genotype column
            genotype = query_row_with_sample['GENOTYPE'].values[0]
            # split on period if this is a double perturbation. Regardless of whether a . is present, genotype will be cast to a list eg ['CNAG_00000'] or ['CNAG_05420', 'CNAG_01438']
            genotype = genotype.split('.')
            # if not wildtype, add to igv_snapshot_dict
            if not genotype[0] == self.wildtype:
                genotype_list.extend(genotype)
                # create bamfile name
                bamfile = sample.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
                # add to igv_snapshot_dict
                igv_snapshot_dict.setdefault(sample, {}).setdefault('gene', []).extend(genotype)
                igv_snapshot_dict[sample]['bam'] = bamfile
                igv_snapshot_dict[sample]['bed'] = None
            # if genotype is equal to wildtype, then store the sample as the wildtype (only one, check if this is right)
            else:
                wt_sample = sample
        # if the wt genotype was found, create entry in the following form {sample_read_counts.tsv: {'gene': [perturbed_gene_1, perturbed_gene_2, ...], 'bam': wt.bam, 'bed': created_bed.bed}
        if wt_sample:
            igv_snapshot_dict[wt_sample]['gene'] = genotype_list
            igv_snapshot_dict[wt_sample]['bam'] = bamfile
            igv_snapshot_dict[wt_sample]['bed'] = None









    def create_igv_region(ineffmut_dict, gene_annot, igv_output_dir, flank, fig_format='png'):
        """
        Create bed files to describe IGV region of interest
        """
        ## get gene dictionary with chromsom, gene coordinates, strand
        if gene_annot.endswith('gtf'):
            gene_annot_dict = utils.parse_gtf(gene_annot)
        elif gene_annot.endswith('gff') or gene_annot.endswith('gff3'):
            gene_annot_dict = utils.parse_gff3(gene_annot)
        else:
            sys.exit("ERROR: The gene annotation format cannot be recognized.")
        ## create gene body region bed file
        for sample in ineffmut_dict.keys():
            igv_bed_filepath = igv_output_dir + sample + '.bed'
            ineffmut_dict[sample]['bed'] = igv_bed_filepath
            writer = open(igv_bed_filepath, 'w')
            for gene in ineffmut_dict[sample]['gene']:
                d = gene_annot_dict[gene]
                writer.write('%s\t%d\t%d\t[%s]%s.%s\n' % \
                             (d['chrm'], d['coords'][0] - flank, d['coords'][1] + flank, sample, gene, fig_format))
        return ineffmut_dict


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

    def index_bams(ineffmut_dict, exp_dir_path):
        """
        Index the bam files in the same directory of bam, and add bam filepath
        """
        # create file in the experiment directory to store the alignment files copied from /lts (these must out of /lts to be accessible to the scheduler)
        experiment_alignment_dir = os.path.join(exp_dir_path, 'alignment_files')
        os.system("mkdir -p {}".format(experiment_alignment_dir))
        # index bams in ineffmut_dict
        for sample in ineffmut_dict.keys():
            # TODO: bam_filepath
            # get the path to the alignment file in /lts
            seq_dir_bam_file = os.path.join('/lts/mblab/Crypto/rnaseq_data',
                                            utils.fileBaseName(sample) + "_sorted_aligned_reads.bam")
            # move the file to the scratch space experiment directory
            os.system("rsync -aHv {} {}".format(seq_dir_bam_file, experiment_alignment_dir))
            # store the path to the bam in scratch
            bam_file = os.path.join(experiment_alignment_dir, os.path.basename(seq_dir_bam_file))
            ineffmut_dict[sample]['bam'] = bam_file
            # index the bam file if it is not already
            if not os.path.isfile(bam_file + '.bai'):
                print('\t Indexing', bam_file)
                out = pysam.index(bam_file)
        return ineffmut_dict