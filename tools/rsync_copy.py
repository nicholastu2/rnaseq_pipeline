#!/usr/bin/env python3
"""
   a script to walk a user through using rsync to move directory and files around
   usage: rsync_copy.py -s <path/to/source> -d <path/to/destination> # formatting doesn't matter much -- the script will ask questions
   to clarify intention and add/remove forward slashes. However, the user SHOULD BE familiar with rsync at least enough to understand
   that the script is asking about whether or not to add or remove forward slashes and why it is important.
"""
import sys
import argparse
from rnaseq_tools import utils
import os


def main(argv):
    """

    :param argv: cmd line arguments
    """
    print(
        '\nWARNING: IF YOU MOVING THE OUTPUT OF ALIGN_COUNTS TO LTS_ALIGN_EXPR, PLEASE USE SCRIPT MOVE_ALIGNMENT_COUNT_FILES.PY\n')
    print('IMPORTANT: read the print statements carefully. Details are important in this one.')
    args = parseArgs(argv)
    print('\nDo you want to copy a directory or individual file(s)? Enter \'d\' or \'f\': ')
    response = input()
    try:
        if response not in ['d', 'f']:
            raise ValueError('MustEnterRecognizedLetter')
    except ValueError:
        print(
            '\nlast chance: only \'d\' or \'f\' are recognized.\nDo you want to copy a directory or individual file(s)? Enter \'d\' or \'f\': ')
        response = input()
        if response not in ['d', 'f']:
            sys.exit('only \'d\' or \'f\' are recognized. Try again from the beginning.')
    if response == 'd':
        source = utils.removeForwardSlash(args.source)
    elif response == 'f':
        print(
            '\nIf this is a single file, enter \'s\'. Else, hit enter. The assumption if you do not enter \'s\' is that\n'
            'you wish to move the contents of a directory. The script will take care of the forward slash formatting: ')
        response = input()
        if response == 's':
            source = args.source
        else:
            source = utils.addForwardSlash(args.source)

    try:
        if not os.path.isdir(args.destination):
            raise FileNotFoundError('DirectoryDoesNotExist')
    except FileNotFoundError:
        print(
            '\nThe Directory you wish to copy the files to Does Not Exist. If you wish to create the directory, enter \'y\'. Else, the script will exit.\n')
        response = input()
        if response == 'y':
            utils.mkdirp(args.destination)
            destination = utils.addForwardSlash(args.destination)
        else:
            sys.exit('Script exiting -- Correct the filepath and try again if you wish.')
    else:
        destination = utils.addForwardSlash(args.destination)

    cmd = 'rsync -aHv %s %s' % (source, destination)
    print('\nexecuting %s\n' % cmd)
    utils.executeSubProcess(cmd)
    print('\nRsync Complete!')


def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="copy with rsync")
    parser.add_argument("-s", "--source", required=True,
                        help='the source of the copy (either a directory, all files in the directory, of a file. Don\'t'
                             'worry about the forward slash and do not add a *. The script will ask you question to figure'
                             'all of this out.')
    parser.add_argument('-d', '--destination', required=True,
                        help='The destination of the copy (a directory)')

    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)
