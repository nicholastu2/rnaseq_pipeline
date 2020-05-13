#!/usr/bin/env python3
"""
   create .bed from gtf with columns chr\tstart\tstop\tfeature
   usage: make_bed.py -a /path/to/.gtf -g CNAG_00102 -b bed_name.bed
"""
import sys
import argparse
import logging
import re


def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # parse cmd line arguments
    args = parseArgs(argv)

    gtf = args.gtf
    gene = args.gene
    bed_name = args.bed_file_name

    file = open(gtf, 'r')
    with open(bed_name, 'w') as bed:
        for line in file:
            if re.search(gene, line):
                x = line.strip().split('\t')
                bed.write('%s\t%s\t%s\t%s\n' % (x[0], x[3], x[4], x[2]))


def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="create a bed file from a gtf and a given gene")
    parser.add_argument("-a", "--gtf", required=True,
                        help='the .gtf for a given organism')
    parser.add_argument('-g', '--gene', required=True,
                        help='the gene you wish to extract from the gtf'),
    parser.add_argument('-b', '--bed_file_name', required=True,
                        help='the name of the bedfile')

    args = parser.parse_args(argv[1:])
    return args


def createLogger(filename, logger_name, logging_conf=None):
    """
    create logger in filemode append and format name-levelname-message with package/module __name__ (best practice from logger tutorial)
    :param filename: name of the file in which to log.
    :param logger_name: __name__ is recommended as a best practice in logger.config eg you can call this like so: createLogger(<your_filename>, __name__)
                     (__name__ is a special variable in python)
    :param logging_conf: path to logging configuration file
    :returns: an instance of the configured logger
    """
    # a config file is passed, load it
    if logging_conf:
        logging.config.fileConfig(logging_conf)  # should include at least what is below
    # if it is not, configure as follows
    else:
        # create log for the year-month-day
        logging.basicConfig(
            filename='%s' % filename,
            filemode='w',
            format='%(name)s-%(levelname)s-%(asctime)s-%(message)s',
            datefmt='%I:%M:%S %p',  # set 'datefmt' to hour-minute-second AM/PM
            level='DEBUG'
        )
    # return an instance of the configured logger
    return logging.getLogger(logger_name)


if __name__ == "__main__":
    main(sys.argv)
