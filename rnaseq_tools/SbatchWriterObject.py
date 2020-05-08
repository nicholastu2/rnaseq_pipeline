
class SbatchWriter:
    def __init__(self, script_path, sbatch_dict, module_list, sample_list, sample_file_extension, command_list):
        self.script_path = script_path
        if sbatch_dict == {}:
            self.sbatch_dict = self.createSbatchDict()
        self.module_list = module_list
        self.sample_list_command = sample_list
        self.sample_file_extension = sample_file_extension
        self.command_list = command_list

    @staticmethod
    def createSbatchDict():
        """
        create a dictionary of sbatch options. Generally, see https://slurm.schedmd.com/sbatch.html
        full names of options are preferred to abbreviations

        :returns: sbatch_dictionary with structure {--sbatch_option: value}
        """
        # create empty dictionary with keys corresponding to sbatch options
        sbatch_options = ['--nodes', '--cpus_per_task', '--mem', '--array', '--chdir', '--output', '--error', '--job-name']
        sbatch_dict = dict.fromkeys(sbatch_options)
        # print message to user
        print('The sbatch options currently available in SbatchWriter.createSbatchDict() are: %s\n'
              'See https://slurm.schedmd.com/sbatch.html for details. Feel free to add more default options'
              'to the script.\n You may also add more to the dictionary after it has been created with this function'
              'and before calling writeScript')
        # prompt user to enter values for each key
        for key in sbatch_dict:
            value = input('enter value for %s: ' %key)
            sbatch_dict[key] = value

        return sbatch_dict

    def writeScript(self):
        with open(self.script_path, 'w') as script_file:
            # write she-bang
            script_file.write('#!/bin/bash\n')
            # write SBATCH options from dict
            for key, value in self.sbatch_dict.items():
                script_file.write('#SBATCH %s=%s\n' %(key, value))
            script_file.write('\n')
            # write modules, if module_list is not empty
            if self.module_list:
                for module in self.module_list:
                    script_file.write('ml %s\n' %module)
            script_file.write('\n')
            if self.sample_list:
                # if a sample_list is passed, create a environmental variable sample_list
                script_file.write("read sample_list < <( sed -n ${{SLURM_ARRAY_TASK_ID}}p {} ); set -e\n\n" %self.sample_list)
                for command in self.command_list:
                    # in the command list, ${{sample}} may be used to write the sample basename as part of a filename. eg:
                    # novoalign -c 8 -o SAM -d {0} -f ${{fastq_file}} 2> {1}/${{sample}}_novoalign.log | samtools view -bS > {1}/${{sample}}_aligned_reads.bam.format(genome_file_index, output_path)
                    script_file.write("sample=${{sample_list##*/}}; sample=${{sample%.%s}}; %s\n" %(self.sample_file_extension,command))
