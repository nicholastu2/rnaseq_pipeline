#!/usr/bin/env python3
"""
   create a lookup file and sbatch script from a user inputted list of runs
   REQUIREMENTS: the runs be in rnaseq_pipeline/align_count_results
   usage: submit_quality_assess_1_batch.py -r 0648 2244 4148
   written by: chase.mateusiak@gmail.com
"""
import sys
import os
import argparse
from rnaseq_tools import utils
from rnaseq_tools.StandardDataObject import StandardData


def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # parse cmd line arguments
    args = parseArgs(argv)
    print('...parsing cmd line arguments')
    query_sheet_path = args.query_sheet
    try:
        if not os.path.isfile(query_sheet_path):
            raise FileNotFoundError('DNE: %s' %query_sheet_path)
    except FileNotFoundError:
        print('The query sheet path is not valid. Check and try again')
    else:
        query_df = utils.readInDataframe(query_sheet_path)

    # store interactive flag
    try:
        interactive_flag = args.interactive
    except AttributeError:
        interactive_flag = False


    run_list = list(query_df.runNumber.unique())

    # create paths from /scratch to the run directory
    sd = StandardData(config_file=args.config_file, interactive=interactive_flag)
    run_path_list = [os.path.join(sd.align_count_results, 'run_'+str(x)+'_samples') for x in run_list]

    # check that paths exist TODO: CHECK CONTENTS OF SUBDIRECTORY FOR COMPLETENESS
    print('...validating paths to run directories')
    validated_run_path_list = validatePaths(sd, run_list, run_path_list)

    # write lookup file of run number paths for the sbatch cmd (see https://htcfdocs.readthedocs.io/en/latest/runningjobs/)
    lookup_filename = 'qual_assess_1_lookup_' + str(sd.year_month_day) + '_' + str(utils.hourMinuteSecond()) + '.txt'
    lookup_output_path = os.path.join(sd.job_scripts, lookup_filename)
    print('...writing lookup file for sbatch script to: %s' %lookup_output_path)
    with open(lookup_output_path, 'w') as file:
        file.write('\n'.join(map(str, validated_run_path_list)))

    # write sbatch script to run qual_assess on all runs in lookup file above
    script = writeSbatchScript(sd, args.user_name, validated_run_path_list, lookup_output_path, query_sheet_path)
    sbatch_filename = 'qual_assess_1_batch_' + str(sd.year_month_day) + '_' + str(utils.hourMinuteSecond() + '.sbatch')
    qual_assess_job_script_path = os.path.join(sd.job_scripts, sbatch_filename)
    print('...writing sbatch script to: %s' %qual_assess_job_script_path)
    with open(qual_assess_job_script_path, "w") as f:
        f.write(script)
    cmd = 'sbatch %s' %qual_assess_job_script_path
    utils.executeSubProcess(cmd)
    print('\nCheck status by cat\'ing the sbatch file above and then cat\'ing the .out file in the sbatch script')

def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="create a lookup file and sbatch script from a user inputted list of runs")
    parser.add_argument("-qs", "--query_sheet", required=True,
                        help="[REQUIRED] A query sheet containing exactly the runs you wish to QA")
    parser.add_argument("-u", "--user_name", required=True,
                        help="[REQUIRED] this MUST BE EXACTLY your htcf user name. you can use the variable $USER to be sure you have it right")
    parser.add_argument('--config_file', default='/see/standard/data/invalid/filepath/set/to/default',
                        help="[OPTIONAL] default is already configured to handle the invalid default path above in StandardDataObject.\n"
                             "Use this flag to replace that config file")
    parser.add_argument('--interactive', action='store_true',
                        help="[OPTIONAL] set this flag (only --interactive, no input necessary) to tell StandardDataObject not\n"
                             "to attempt to look in /lts if on a compute node on the cluster")

    args = parser.parse_args(argv[1:])
    return args

def validatePaths(sd, run_list, run_path_list):
    """
        check that run files exist where expected, in rnaseq_pipeline/align_count_results
        :param sd: a standard data object
        :param run_path_list: a list of full paths from / to run_#####_samples created in main script
        :returns: a validated list
    """
    # check that the run directories all exist
    for i in range(len(run_path_list)):
        path = run_path_list[i]
        try:
            if not os.path.isdir(path):
                try:
                    run_num_with_leading_zero = sd._run_numbers_with_zeros[run_list[i]]
                except KeyError:
                    sd.logger.info('%s not in leading zero run number list')
                    raise NotADirectoryError('DoesNotExist: %s' %path)
                else:
                    path_with_zero = os.path.join(sd.align_count_results, 'run_'+str(run_num_with_leading_zero)+'_samples')
                    if not os.path.isdir(path_with_zero):
                        raise NotADirectoryError('DoesNotExist: %s' %path)
                    else:
                        run_path_list[i] = path_with_zero
        except NotADirectoryError:
            msg='Cant find the run directory. The runnumbers must be subdirectories like so: rnaseq_pipeline/align_count_results/run_####_samples'
            sd.logger.info(msg)
            print(msg)

        return run_path_list

def writeSbatchScript(sd, user_name, validated_run_path_list, lookup_output_path, query_path):
    """
        write sbatch script to qual assess runs
        :param sd: a standard data object
        :param validated_run_path_list: list of run paths, for which the paths have been verified as existing
        :param lookup_output_path: path to the lookup file (see writeLookupFile())
        :param query_path: path to query
        :returns: a formatted sbatch script stored as a string
    """
    script= ("#!/bin/bash\n"
             "#SBATCH -N 1\n"
             "#SBATCH --cpus-per-task=8\n"
             "#SBATCH --mem=12G\n"
             "#SBATCH --array=1-{0}%{1}\n".format(len(validated_run_path_list), min(len(validated_run_path_list), 50)))
    script = script + "#SBATCH -D /scratch/mblab/%s/rnaseq_pipeline\n" \
                      "#SBATCH -o sbatch_log/qual_assess_%s.out\n" \
                      "#SBATCH -J qual_assess_1_batch\n\n" %(user_name, str(sd.year_month_day) + '_' + str(utils.hourMinuteSecond()))
    script = script + 'ml rnaseq_pipeline\n\n'
    script = script + 'read run_path < <( sed -n ${SLURM_ARRAY_TASK_ID}p %s )\n\n' %lookup_output_path
    script = script + 'quality_assess_1.py -ac ${run_path} -qs %s --interactive\n' %query_path

    return script


if __name__ == "__main__":
    main(sys.argv)