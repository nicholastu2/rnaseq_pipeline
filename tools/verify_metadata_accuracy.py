#!/usr/bin/env python

# name: verify_metadata_accuracy
# purpose: check metadata tables for expected entries
# input: an individual sheet, a directory of sheets, or the top most directory
# output: Warnings and prompts to the user to either ignore or change entries
# written by: chase mateusiak(chase.mateusiak@gmail.com)
# date included in rnaseq_pipe: 1/21/2020
# See the end of the script for a description of the environment used to write/test

#import queryDB
import pandas as pd
import sys
import argparse
import re
import os
from queryDB import *

unique_key_columns = {"fastqFiles": ['libraryDate', 'libraryPreparer', 'librarySampleNumber'],
        "library":['libraryDate', 'libraryPreparer', 'librarySampleNumber','s2cDNADate', 's2cDNAPreparer', 's2cDNASampleNumber'],
        "s2cDNASample": ['s2cDNADate', 's2cDNAPreparer', 's2cDNASampleNumber','s1cDNADate', 's1cDNAPreparer', 's1cDNASampleNumber'],
        "s1cDNASample": ['s1cDNADate', 's1cDNAPreparer', 's1cDNASampleNumber', 'rnaDate', 'rnaPreparer', 'rnaSampleNumber'],
        "rnaSample": ['rnaDate', 'rnaPreparer', 'rnaSampleNumber','harvestDate', 'harvester', 'biosampleNumber'],
        "bioSample": ['harvestDate', 'harvester', 'biosampleNumber']}

def main(argv):
    args = parseArgs(argv)

    if args.database_subdir:
        datadir_dict = getFilePaths(args.database, args.database_subdir)
    else:
        datadir_dict = getFilePaths(args.database)

    concat_dict = {}

    for key in datadir_dict:
        concat_dict[key] = concatMetadata(datadir_dict[key])
        cols = '\t'.join(concat_dict[key].columns)
        cols = cols.lower()
        if "unnamed" in cols:
            print("\nThere is an unnamed column in one of the sheets in subdirectory {}. This must be found and fixed before proceeding.".format(key))
            quit()

    for key in concat_dict:
        uniqueKeys(unique_key_columns[key], concat_dict[key])

    print('\n::verify_metadata_complete::\nPlease fix any issues found before pushing to remote')

def parseArgs(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', required=True,
                        help='Suggested usage: database-files path to metadata base')
    parser.add_argument('-k', '--database_subdir', required=False,
                        help='if you do not wish to verify all metadata subdirectories, you may list the ones you would like to check. Default is to check all metadata subdirs in database')

    return parser.parse_args(argv[1:])

def checkCSV(file):
    # test whether a given file is a .csv or .xlsx
    if re.search('\.csv', file):
        return True
    else:
        return False

def getKeys(datadir_keys, concat_dict):
    # create list of shared columns between successive pairs of keys in datadir_keys i.e. the columns which are shared between the sheets in fastqFiles and Library
    # Args: a list of subdirectories in datadir and a dictionary of concatenated sheets from within each one of those directories
    # PLEASE NOTE: DICTIONARIES ARE ORDERED IN PYTHON 3.6+. KEY/DIRECTORIES NEED TO BE ADDED TO CONCAT DICT IN ORDER YOU WISH TO MERGE
    # Returns: a list of key columns
    key_cols = []
    num_keys = len(datadir_keys)
    for i in range(num_keys-1):
        key_cols.append(concat_dict[datadir_keys[i]].keys().intersection(concat_dict[datadir_keys[i+1]].keys()))
    return key_cols

def uniqueKeys(key, df):
    # test whether the key columns in a given sheet are unique
    # Args: the keys (passed as list) and the path to the dataframe
    # Return: none
    # to std_out: info on uniqueness of key
    if not isinstance(df, pd.DataFrame) and os.path.isfile(df):
        if checkCSV(df):
            sheet = pd.read_csv(df)
        else:
            sheet = pd.read_excel(df)
    else:
        sheet = df
    # make a tuple of the key columns and store them as pandas series
    key_tuples = sheet[key].apply(tuple, axis=1)
    num_keys = key_tuples.size
    num_unique_keys = key_tuples.unique().size
    print("\nThe number of unique keys is {}. The number of rows is {}. If these are equal, the keys are unique.".format(num_keys, num_unique_keys))
    if not num_keys == num_unique_keys:
        print("\nThe following indicies are not unique:\n\t{}".format(sheet[sheet[key].duplicated()]))
        print("\nDo you want to continue? Enter y or n: ")
        user_response = input()
        if user_response == 'n':
            quit()

def coerceAllCols(sheet):
    # converts columns specific in functions below to floats and datetime
    sheet = floatColCoerce(sheet)

    sheet = datetimeColCoerce(sheet)

    # forcing all to lowercase is not working in the search
    #sheet = strColCoerce(sheet)

    return sheet

def colCoerce(cols, sheet, dtype):
    # general function to coerce column to certain data type
    # Args: a list of columns to coerce (can be a subset of columns of a data frame. must be a list), a dataframe,
    #       and a datatype (see pandas .astype documentation for which datatypes are possible. didn't work for datetime, for example)
    # Return: the dataframe with the specified columns coerced
    coerce_cols = cols

    for col in sheet.columns[sheet.columns.isin(coerce_cols)]:
        sheet[col] = sheet[col].astype(dtype)
    return sheet

def floatColCoerce(sheet):
    # coerce columns specified below to int
    int_columns = ['librarySampleNumber', 'runNumber', 'tapestationConc', 'readsObtained', 'rnaSampleNumber',
                   'inductionDelay', 'replicate', 'timePoint']

    sheet = colCoerce(int_columns, sheet, 'float')

    return sheet

def datetimeColCoerce(sheet):
    # coerce columns specified below to datetime
    datetime_columns = ['libraryDate', 's2cDNADate', 's1cDNADate', 'rnaDate', 'harvestDate']

    for col in datetime_columns:
        sheet[col] = pd.to_datetime(sheet[col])

    return sheet

def strColCoerce(sheet):
    # coerce any column NOT in int_col or datetime_col to a string and make it lower case

    int_columns = ['librarySampleNumber', 'tapestationConc', 'readsObtained', 'rnaSampleNumber','inductionDelay', 'replicate', 'timePoint']
    datetime_columns = ['libraryDate', 's2cDNADate', 's1cDNADate', 'rnaDate', 'harvestDate']

    typed_cols = int_columns + datetime_columns

    for col in sheet.columns:
        if col not in typed_cols:
            sheet[col] = sheet[col].astype(str).apply(lambda x: x.lower())

    return sheet

if __name__ == '__main__':
	main(sys.argv)
