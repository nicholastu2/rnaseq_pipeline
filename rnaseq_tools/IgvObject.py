"""
tools to create browser shots

makeIgvSnapshotDict create a dict with the following format
igv_snapshot_dict = {'wildtype.fastq.gz': {'gene': [gene_1, gene_2, gene_3...], 'bam': file.bam, 'bed': file.bed},
                     'perturbed_sample.fastq.gz': {'gene': [gene_1, gene_2], 'bam': file.bam, 'bed': file.bed},
                     'perturbed_sample.fastq.gz': {'gene': [gene_1], 'bam': file.bam, 'bed': file.bed}
                    }
and ensures that all alignment files (.bam) and their index companions (.bam.bai) are in the experiment directory, which
is either supplied or created in rnaseq_tmp

usage: igv = IgvObject(organism = <your organism>, query_sheet_path = <path to query>,,
                       igv_output_dir = output directory for igv snapshots,
                       scratch_alignment_source = <path to experiment dir> # if this is not present, bed/bat files will be put in rnaseq_tmp/datetime_igv_files)

"""

from rnaseq_tools import utils
from rnaseq_tools.OrganismDataObject import OrganismData
from rnaseq_tools import annotation_tools
from rnaseq_tools.DatabaseObject import DatabaseObject
import sys
import os
import time

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
        self._igv_attributes = ['query_sheet_path', 'igv_output_dir']
        # initialize Standard data with the extended _attributes
        super(IgvObject, self).__init__(self._igv_attributes, **kwargs)
        # initialize list to store bamfiles that need to be indexed (must be done by batch script)
        self.bam_file_to_index_list = []
        # create logger for IgvObject
        self.logger = utils.createLogger(self.log_file_path, __name__, 'DEBUG')
        try:
            self.query_sheet_path = kwargs['query_sheet_path']
            self.control_sheet_path = kwargs['control_sheet_path']
        except KeyError:
            self.logger.debug('query sheet path and/or control sheet path not passed in constructor')
        else:
            self.sample_df = utils.readInDataframe(self.query_sheet_path)
            self.control_sample_df = self.createControlSampleDict()

        # get gene dictionary with chromsome, gene coordinates, strand
        if self.annotation_file.endswith('gtf'):
            self.annotation_dict = annotation_tools.parseGtf(self.annotation_file)
        elif self.annotation_file.endswith('gff') or self.annotation_file.endswith('gff3'):
            self.annotation_dict = annotation_tools.parseGff3(self.annotation_file)
        else:
            sys.exit("ERROR: The gene annotation format cannot be recognized.")  # TODO: clean up preceeding blocks -- move parseGFF to OrganismData


    def createControlSampleDict(self): # TODO: put this in KN99 genome files
        """
            from a query sheet filtered for a single wildtype sample per treatment_timepoint ( or some other set of columns )
            create a dictionary with structure:
                {treatment_timepoint: bam_fullpath} eg {30C.CO2_90: sample1_....bam, 37C.CO2: sample2_....bam, ... }
            note: the expected source of the wildtype files is in self.align_count_results/run_####_samples/align/bamfilename.bam
        """
        wildtype_dict = {}

        for index, row in self.wildtype_df.iterrows():
            fastq_filename = row.fastqFileName
            bam_filename = utils.convertFastqFilename(fastq_filename, 'bam')
            run_number = row.runNumber # TODO: this is copied code - either put in function or use utils.extractFromQuerySheet
            try:
                run_num_with_zero = self._run_numbers_with_zeros[run_number]
            except KeyError:
                pass
            else:
                run_number = run_num_with_zero
            treatment = utils.extractInfoFromQuerySheet(self.wildtype_df, fastq_filename, 'treatment')
            timepoint = utils.extractInfoFromQuerySheet(self.wildtype_df, fastq_filename, 'timePoint')
            bam_file_fullpath = os.path.join(self.align_count_results, 'run_%s_samples/align'%run_number, bam_filename)
            try:
                if not os.path.exists(bam_file_fullpath):
                    raise FileNotFoundError('control sample %s not found' %bam_file_fullpath)
            except FileNotFoundError:
                error_msg = 'control sample not in %s. %s DNE' %(self.align_count_results, bam_file_fullpath)
                self.logger.critical(error_msg)
                print(error_msg)
            wildtype_dict.setdefault("%s_%s" %(treatment, timepoint), bam_file_fullpath)

        return wildtype_dict

    def makeIgvSnapshotDict(self, control_sample_flag=True, marker_list=['CNAG_G418', 'CNAG_NAT']):
        """ TODO: TEST THAT .BAI FILE IS PRESENT
        TODO: RIGHT NOW SPECIFICALLY SET UP FOR 90MINUTEINDUCTION CRYPTO (because of control and marker)
            from list of samples, create dictionary with structure
                {fastq_filename1: {'gene': [gene_1, gene_2,...], 'bed': /path/to/bed, 'bam': /path/to/bam, fastq_filename2: ...}
        """
        setattr(self, 'igv_snapshot_dict', {})
        igv_sample_dict = {}

        for index, row in self.sample_df.iterrows():
            # extract relevant info from query row
            sample_name = utils.pathBaseName(row.fastqFileName)
            bam_file = utils.convertFastqFilename(row.fastqFileName, 'bam')
            treatment = str(row.treatment)
            timepoint = str(row.timePoint)
            treatment_timepoint = "%s_%s" %(treatment,timepoint) #NOTE: this is setup specifically for KN99 -- needs to be generalized
            run_number = row.runNumber
            try:
                run_num_with_zero = self._run_numbers_with_zeros[run_number]
            except KeyError:
                pass
            else:
                run_number = run_num_with_zero
            bamfile_fullpath = os.path.join(self.align_count_results, 'run_%s_samples/align' %run_number, bam_file)
            try:
                if not os.path.exists(bamfile_fullpath):
                    raise FileNotFoundError('%s DNE' %bamfile_fullpath)
            except FileNotFoundError:
                error_msg = 'bamfile does not exist at %s' %bamfile_fullpath
                self.logger.critical(error_msg)
                print(error_msg)
            genotype = row.genotype
            # split on period if this is a double perturbation. Regardless of whether a '.' is present,
            # genotype will be cast to a list eg ['CNAG_00000'] or ['CNAG_05420', 'CNAG_01438']
            genotype_list = genotype.split('.')
            # if genotype (possibly two if double KO, hence loop) ends with _over, remove _over
            genotype_list = [x.replace('_over', '') for x in genotype_list]
            # replace any CNAG with CKF44 -- for crypto, CNAG was h99 designation and continues to be used in metadata. CKF44 is kn99 prefix in ncbi and is the prefix used in the genome, annotation, etc
            genotype_list = [x.replace('CNAG', 'CKF44') for x in genotype_list]

            # add sample to igv_sample_dict
            igv_sample_dict.setdefault(sample_name, {})
            if control_sample_flag:
                igv_sample_dict[sample_name]['control_sample'] = self.getControlSample(treatment_timepoint)
            # add locus_dict
            igv_sample_dict[sample_name]['locus_dict'] = self.createLocusDict[genotype_list]
            # add marker_dict --> this will be pretty standard, maybe put in attributes somewhere eventually
            if marker_list:
                igv_sample_dict[sample_name]['marker_locus_dict'] = self.createLocusDict[marker_list]

    def getControlSample(self, key):
        """ TODO: makeIGvDict needs to be generalized so that it doesn't always extract treatment_timepoint and put that into getControlSample (eg, if there is a different set of columns taht describe a control samlpe)
            in preparation for a more general code, this returns the corresponding control sample given a key
            :params treatment: treatment condition
            :params timepoint: timepoint condition
            :returns: bamfilepath, expected source is self.align_count_results
        """
        return self.wildtype_dict[key]

    def createLocusDict(self, genotype_list):
        """
            given a genotype list, return a dictionary like so:
               {genotype_list[0]: [chrm, start, stop], genotype_list[1]: [chrm, start, stop]
        """
        locus_dict = {}

        for locus in genotype_list:
            locus_info = self.annotation_dict[locus]
            locus_dict.setdefault(locus, []).extend([locus_info['chrm'], locus_info['coords'][0], locus_info['coords'][1]])

        return locus_dict

    def createBedFile(self, gene_list, fastq_filename, description, flanking_region=500):
        """ TODO currently file is _read_count.bed -- git rid of the _read_count
            Create bed files to describe IGV region of interest
            :param gene_list: list of genes to include in the bed file
            :param fastq_filename: fastq filename (sample identifier). used to name the bedfile and name the track in the igv shot
            :param description: either 'perturbed' or 'wildtype', or some other description of the track
            :param flanking_region: how far up/down stream from the gene to capture in the snapshot
            :returns: bedfile path (will be in rnaseq_tmp/todays_date
        """

        # create day specific directory in rnaseq_tmp
        rnaseq_tmp_bedfile_dirpath = os.path.join(self.rnaseq_tmp, self.year_month_day)
        utils.mkdirp(rnaseq_tmp_bedfile_dirpath)

        # sample is the fastq file basename
        sample = utils.pathBaseName(fastq_filename)

        # create gene body region bed file
        igv_bed_filepath = os.path.join(rnaseq_tmp_bedfile_dirpath, utils.pathBaseName(fastq_filename) + '.bed')

        bed_lines_list = []
        for gene in gene_list:
            gene_parsed_annotation_dict = self.annotation_dict[gene]
            sample_description = '[%s]%s.png' %(sample, gene)
            bed_lines_list.append('%s\t%s\t%s\t%s\t%s\t%s\n' %(gene_parsed_annotation_dict['chrm'],
                                                                max(gene_parsed_annotation_dict['coords'][0] - flanking_region, 1),
                                                                gene_parsed_annotation_dict['coords'][1] + flanking_region,
                                                                sample_description, ".", gene_parsed_annotation_dict['strand']))

        with open(igv_bed_filepath, 'w') as file:
            file.write(''.join(bed_lines_list))

        return igv_bed_filepath

    @staticmethod
    def writeIgvBatchScript(output_dir, batchfilename, genome, bam_file_list, locus_dict):
        """
        create an IGV batch script (see example and site below)
        :param output_dir: output directory for bathfile
        :param batchfilename: name of batchfile
        :param genome: .genome file created by the igv browser
        :param bam_file_list: list of alignment files to load (eg, to compare a control to a perturbed sample, there would be two items in the list)
                              A singular bam file should still be passed as a list
        :param locus_dict: structure of dict:
            {locus1: {chrm: chr1, loc: [start, stop],
             locus2: ...}
        :returns: batchfile path

            The arguments (unique_batchfilename, genome, [control.bam, perturbed.bam, {locus1: {chrm: chr1, loc: [1,100], locus2: {chrm: chr10, loc: [3,300]}}}])
             will create the following batchscript:
        http://software.broadinstitute.org/software/igv/book/export/html/189
            example IGV batch script:
                new                                         # batchscript keyword new (new snapshot)
                snapshotDirectory IGV_Snapshots             # output directory
                maxPanelHeight 500                          # maximum height of igv browser viewer
                genome <input_genome>                       # igv .genome file
                load control.bam                            # load alignment files in bam_file_list
                load perturbed.bam
                goto chr1:1-100                             # load region of interest
                snapshot <batchfilename>_locus1 .png        # saves a snapshot of the IGV window to an image file
                goto chr10:3-300                            # repeat at another locus
                snapshot <batchfilename>_locus2.png
                exit                                        # close igv session

                # Note: in this case, there will be two tracks, the first being the control, the next being the perturbed sample, one on top of the other
                        If only one locus is passed, there will be only one track (plus the annotation track) in view
        """

    def writeIgvJobScript(self, email = None, fig_format='png'):
        """
        Write sbatch job script to make IGV snapshot
           Calls helper function -- write igv_lookup to write lookup file
        :param email: option email address to add to sbatch
        :param fig_format: default figure format is png
        """
        lookup_file_path = self.writeLookupFile()

        # set igv_output_dir if not entered on cmd line
        if not hasattr(self, 'igv_output_dir'):
            setattr(self, 'igv_output_dir', os.path.join(self.scratch_alignment_source, 'igv_output'))
        num_samples = len(self.igv_snapshot_dict.keys())
        job = '#!/bin/bash\n' \
              '#SBATCH -N 1\n' \
              '#SBATCH --array=1-{0}%20\n' \
              '#SBATCH --mem=10G\n' \
              '#SBATCH -o {1}/igv_snapshot_%A_{2}.out\n' \
              '#SBATCH -e {1}/igv_snapshot_%A_{2}.err\n' \
              '#SBATCH -J igv_snapshot_{2}\n'.format(num_samples, self.sbatch_log, self.year_month_day + '_' + utils.hourMinuteSecond())
        if email:
            job += '#SBATCH --mail-type=END,FAIL\n' \
                   '#SBATCH --mail-user=%s\n' % self.email
        job += '\nml java\n' \
               'ml rnaseq_pipeline\n\n' \
               'read bam_file bed_file igv_genome < <(sed -n ${SLURM_ARRAY_TASK_ID}p %s )\n\n' \
               'make_IGV_snapshots.py $bam_file -bin /opt/apps/igv/2.4.7/igv.jar -nf4 -r $bed_file -g $igv_genome -fig_format %s -o %s\n'\
                                                                         %(lookup_file_path, fig_format, self.igv_output_dir)

        igv_job_script = os.path.join(self.today_job_dir, '%s_%s.sbatch' %(self.year_month_day, utils.hourMinuteSecond()))
        with open(igv_job_script, 'w') as file:
            file.write('%s' % job)
        # TODO add a log statement, think about redirecting these job scripts to the job_script directory
        print('the job scripts is at %s\n'
              '%s has been set as the igv_output_directory' % (igv_job_script, self.igv_output_dir))  # check this

    def writeLookupFile(self, lookup_filename = None):
        """
            write a lookup file for igv snapshot sbatch script
            :param lookup_filename: Default None, in which case the name will be lookup_year_month_day_hour_minute_second
            :returns: path to lookup file
        """
        try:
            igv_genome = self.igv_genome
        except AttributeError:
            igv_genome = self.igv_stranded_genome

        if lookup_filename is None:
            lookup_filename = 'lookup_%s_%s.txt' %(self.year_month_day, utils.hourMinuteSecond())

        lookup_file_list = []
        for sample in self.igv_snapshot_dict.keys():
            bam_file = self.igv_snapshot_dict[sample]['bam']
            bed_file = self.igv_snapshot_dict[sample]['bed']
            lookup_file_list.append('%s\t%s\t%s\n' %(bam_file, bed_file, igv_genome))

        setattr(self, 'today_job_dir', os.path.join(self.job_scripts, self.year_month_day))
        utils.mkdirp(self.today_job_dir)

        lookup_file_path = os.path.join(self.today_job_dir, lookup_filename)
        with open(lookup_file_path, 'w') as file:
            file.write(''.join(x for x in lookup_file_list))

        return lookup_file_path

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
