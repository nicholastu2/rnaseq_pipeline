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

def uniqueKeys(key, sheet):
    key = list(key)
    # make a tuple of the key columns and store them as pandas series
    key_tuples = sheet[key].apply(tuple, axis=1)
    print(len(key_tuples), len(sheet.index))

    if not key_tuples.unique().size == len(sheet.index):
        print("\nThe following indicies are not unique:\n\t{}".format(key_tuples.unique))
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

uniqueKeys(['libraryDate','libraryPreparer','librarySampleNumber',
            's2cDNADate', 's2cDNAPreparer', 's2cDNASampleNumber'], pd.read_csv('~/Desktop/first_merge.csv'))