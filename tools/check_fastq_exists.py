#!/usr/bin/env python

import pandas as pd
#from rnaseq_tools.StandardDataObject import StandardData
import os
import argparse
import sys


# TODO: XCODE, WGET, SSHFS
#       FIGURE OUT PERMISSIONS FOR PYCHARM (MAYBE NOT NECESSARY?)
#       get the config_file and interactive arguments working in parseArgs (use the debugger and breakpoints to work through this)
#

def main(argv):
    args = parseArgs(argv)
    path = args.query_sheet
    df = pd.read_csv(path)

    #sd = StandardData(config_file='/Users/tomtu/Documents/brentlab/rnaseq_config.ini')
    #TODO make command line arguments to pass in the path and interactive=True

    def extractRunNumber(sd, run_num):
        try:
            return str(sd._run_numbers_with_zeros[int(run_num)])
        except KeyError:
            return str(int(run_num))

    for index, row in df.iterrows():
        fastq = row.fastqFileName
        #run_num = extractRunNumber(sd, float(row.runNumber))
        # fastq_path = os.path.join(sd.lts_sequence, "run_"+run_num+"_samples", fastq)
        #TODO look up if str works
        if not isinstance(fastq, str):
            print(fastq)
    print("all there!")

def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="For a given alignment/count job, create a nextflow config file")
    parser.add_argument("-qs", "--query_sheet", required=True,
                        help='[REQUIRED] A .csv subset of metadata database describing the set of files you wish to process.\n'
                             'These will be grouped by run number')#TODO add config file argument

    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)

