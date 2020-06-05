import os

class SbatchWriter:
    def __init__(self, **kwargs):
        try:
            self.script_path = kwargs['script_path']
        except KeyError:
            pass
        try:
            self.sbatch_dict = kwargs['sbatch_dict']
        except KeyError:
            pass
        try:
            self.module_list = kwargs['module_list']
        except KeyError:
            pass
        try:
            self.sample_list = kwargs['sample_list']
        except KeyError:
            pass
        try:
            self.sample_file_path = kwargs['sample_file_path']
        except KeyError:
            pass
        try:
            self.command_list = kwargs['command_list']
        except KeyError:
            pass

    # TODO: test and then incorporate SbatchWriterObject
    @staticmethod
    def writeAlignCountJobScript(job_file, output_path, fastq_list_file, num_fastqs, genome_index_file,
                                 genome_annotation_file, feature_type, strandness, align_only):
        """
        Write slurm job script to job_file (which is $PWD/job_scripts
        :param job_file: path to $PWD/job_scripts (see main method)
        :param output_path: path to output_dir (see main method)
        :param fastq_list_file: path to job_scripts/<run_number>_fastq_list.txt
        :param num_fastqs: number of lines in fastq_list_file
        :param genome_index_file: path to index file created from genome of given organism by novoalign
        :param genome_annotation_file: path to genome annotation file, either .gff or .gtf
        :param feature_type: feature type extracted based on FEATURE_TYPE_DICT (see top of script)
        :param strandness: cmd line input from user regarding whether library prep is stranded
        :param align_only: boolean flag allowing cmd line input for alignment only, no htseq
        :returns: None
        """

        with open(job_file, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH -N 1\n")
            f.write("#SBATCH --cpus-per-task=8\n")
            f.write("#SBATCH --mem=12G\n")
            f.write("#SBATCH --array=1-{0}%{1}\n".format(num_fastqs, min(num_fastqs, 50)))
            f.write("#SBATCH -D ./\n")
            f.write("#SBATCH -o sbatch_log/mblab_rnaseq_%A_%a.out\n")
            f.write("#SBATCH -e sbatch_log/mblab_rnaseq_%A_%a.err\n")
            f.write("#SBATCH -J mblab_rnaseq\n\n")

            f.write("ml novoalign/3.07.00\n")
            f.write("ml samtools/1.6\n")
            f.write("ml htseq/0.9.1\n")
            f.write("read fastq_file < <( sed -n ${{SLURM_ARRAY_TASK_ID}}p {} ); set -e\n\n".format(fastq_list_file))
            f.write("mkdir -p {}\n".format(output_path))

            f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; novoalign -r All -c 8 -o SAM -d {0} "  # NOTE: ADDED -r all
                    "-f ${{fastq_file}} 2> {1}/${{sample}}_novoalign.log | samtools view -bS > "
                    "{1}/${{sample}}_aligned_reads.bam\n".format(genome_index_file, output_path))

            f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; novosort --threads 8 "
                    "{0}/${{sample}}_aligned_reads.bam > {0}/${{sample}}_sorted_aligned_reads.bam 2> "
                    "{0}/${{sample}}_novosort.log\n".format(output_path))

            if not align_only:
                if feature_type == 'gene':  # this is a messy way of saying "if gff, specify -i ID. TODO: This needs to be cleaned up
                    f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; "
                            "htseq-count -f bam -i ID -s {1} -t {2} {0}/${{sample}}_sorted_aligned_reads.bam "
                            "{3} > {0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path,
                                                                                                         strandness,
                                                                                                         feature_type,
                                                                                                         genome_annotation_file))
                else:  # else (it is a gtf) do not specify ID
                    f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; htseq-count -f bam -s {1} "
                            "-t {2} {0}/${{sample}}_sorted_aligned_reads.bam {3} > "
                            "{0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path,
                                                                                                   strandness,
                                                                                                   feature_type,
                                                                                                   genome_annotation_file))

    @staticmethod
    def writeCoverageSbatchScript(job_file, output_dir_path, run_num_directory, perturbed_sorted_alignment_list):
        """

        """
        with open(job_file, 'w') as sbatch_file:
            sbatch_file.write("#!/bin/bash\n")
            sbatch_file.write("#SBATCH --mem=5G\n")
            sbatch_file.write("#SBATCH -D ./\n")
            sbatch_file.write("#SBATCH -o sbatch_log/coverage_calculation_%A_%a.out\n")
            sbatch_file.write("#SBATCH -e sbatch_log/coverage_calculation_%A_%a.err\n")
            sbatch_file.write("#SBATCH -J coverage_calculation\n\n")
            sbatch_file.write("ml bedtools\n\n")
            for sorted_alignment_file in perturbed_sorted_alignment_list:
                sorted_alignment_path = os.path.join(run_num_directory, sorted_alignment_file)
                coverage_file_name = sorted_alignment_file.replace('_sorted_aligned_reads.bam ', 'coverage.tsv')
                coverage_file_path = os.path.join(output_dir_path, coverage_file_name)
                sbatch_file.write('bedtools genomecov -ibam %s -bga > %s\n' % (sorted_alignment_path, coverage_file_path))

    @staticmethod
    def createSbatchDict():
        """
        create a dictionary of sbatch options. Generally, see https://slurm.schedmd.com/sbatch.html
        full names of options are preferred to abbreviations

        :returns: sbatch_dictionary with structure {--sbatch_option: value}
        """
        # create empty dictionary with keys corresponding to sbatch options
        sbatch_options = ['--nodes', '--cpus_per_task', '--mem', '--array', '--chdir', '--output', '--error',
                          '--job-name']
        sbatch_dict = dict.fromkeys(sbatch_options)
        # print message to user
        print('The sbatch options currently available in SbatchWriter.createSbatchDict() are: %s\n'
              'See https://slurm.schedmd.com/sbatch.html for details. Feel free to add more default options'
              'to the script.\n You may also add more to the dictionary after it has been created with this function'
              'and before calling writeScript')
        # prompt user to enter values for each key
        for key in sbatch_dict:
            value = input('enter value for %s: ' % key)
            sbatch_dict[key] = value

        return sbatch_dict

    def writeSbatchHeader(self):
        with open(self.script_path, 'w') as script_file:
            # write she-bang
            script_file.write('#!/bin/bash\n')
            # write SBATCH options from dict
            for key, value in self.sbatch_dict.items():
                script_file.write('#SBATCH %s=%s\n' % (key, value))
            script_file.write('\n')
            try:
                for module in self.module_list:
                    script_file.write('ml %s\n' % module)
                script_file.write('\n')
            except AttributeError:
                pass
