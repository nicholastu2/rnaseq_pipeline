#!/usr/bin/env python

# name: populate_fastqfileName
# purpose: populate field: fastqFileName of fastqFile directory metadata sheets
# input: topmost database directory of rnaseq metadata, directory containing fastq (or fastq.gz etc) files
# output: fastqFile metadata sheets populated with fastqFileNames (where the fastq file exists)
# written by: chase mateusiak(chase.mateusiak@gmail.com)
# date included in rnaseq_pipe: 1/20/2020
# See the end of the script for a description of the environment used to write/test

# TODO: ignore case (possibly this function lives in different script. maybe the accuracy script prior to get push database)

# import necessary functions from queryDB in ./tools
from queryDB import getFilePaths, createDB, checkCSV
import pandas as pd
import fnmatch
import os
import re

def checkCSV(file):
    # test whether a given file is a .csv or .xlsx
    if re.search('\.csv', file):
        return True
    else:
        return False

def createFastqMetadata(metadata_dir, datadir_keys):
    # creates dataframe of concatenated metadata sheets in fastqFiles and library
    # Args: top level directory of metadata; subdirectories to search for data sheets. In this case, fastqFiles and library
    # Returns: a datasheet with information necessary to parsing the sequence directory

    datadir_dict = getFilePaths(metadata_dir, datadir_keys)
    concat_df = createDB(datadir_dict, datadir_keys, drop_fastq_na=False)
    #concat_df.to_csv("~/Desktop/fastq_parsing.csv")
    return concat_df

def fastqInput(fastqFile_metadata, metadata_df, sequence_dir):
    #TODO: Currently, writes out sheets with "nan" for unfilled values. Should this be changed?

    # populates fastqFileName field of fastqFile sheets
    # Args: fastqFiles is EITHER the topmost metadata directory OR the path to a single sheet.

    #       To update filepaths for MORE THAN ONE, BUT LESS THAN ALL fastqFile sheets, create a
    #       separate directory (name/location doesn't matter) with subdirectory 'fastqFiles'
    #       place the batch of sheets you'd like to update in the subdir

    # Return: None. The given fastqFiles sheets are updated with fastqFileName IN PLACE

    # get list of filepaths (getFilePaths returns as dictionary of lists. extract list)
    if os.path.isdir(fastqFile_metadata):
        fastqFiles_paths = getFilePaths(fastqFile_metadata, ['fastqFiles'])
        fastqFiles_paths = fastqFiles_paths['fastqFiles'][0]
    else:
        fastqFiles_paths = fastqFile_metadata

    # if list, iterate through the sheets in the fastqFile_Paths
    if isinstance(fastqFiles_paths, list):
        for sheet_path in fastqFiles_paths:
            if checkCSV(sheet_path):
                sheet_df = pd.read_csv(sheet_path)
            else:
                sheet_df = pd.read_excel(sheet_path, index_col=None)
            # if not list
            sheet_df = populateFastqFileName(sheet_df, metadata_df, sequence_dir)
            sheet_df.to_excel(sheet_path, na_rep='', index=False)

    else:
        # read in sheet
        if checkCSV(fastqFiles_paths):
            sheet_df = pd.read_csv(fastqFiles_paths)
        else:
            sheet_df = pd.read_excel(fastqFiles_paths)
        sheet_df = populateFastqFileName(sheet_df, metadata_df, sequence_dir)
        sheet_df.to_excel(fastqFiles_paths, na_rep='', index=False)

def populateFastqFileName(sheet_df, metadata_df, sequence_dir):
    #TODO: Do not enter if less than certain number of reads? Don't overwrite without user prompt; if a fastq filename exists, check if it matches
        # THIS NEEDS TO BE DONE AS PART OF A "VERIFY_STRUCTURAL_ACCURACY" SCRIPT

    # coerce sheet to all strings
    sheet_df = sheet_df.astype('str')

    # iterate through the unique run numbers of a given sheet (though typically it seems there is only one unique run number)
    for runNum in sheet_df.runNumber.unique():

        # store the expected name of the raw fastq data directory as run_[runNum]_samples
        fastq_dir = 'run_' + str(runNum) + '_samples'

        # test whether the raw fastq file exists. If so, return unaltered sheet
        fastq_dir_full_path = os.path.join(sequence_dir, fastq_dir)
        if not os.path.exists(fastq_dir_full_path):
            return sheet_df
        else:
            # iterate through the rows of the datasheet by unique run number
            for index, row in sheet_df[sheet_df.runNumber == runNum].iterrows():
                # variables to identify a row uniquely
                lib_date = row['libraryDate']
                lib_prepper = row['libraryPreparer']
                sample_num = row['librarySampleNumber']
                purpose = row['purpose']

                # create query formula from the variables above
                query_formula = f'libraryDate == "{lib_date}" & libraryPreparer == "{lib_prepper}" & runNumber == "{runNum}" & librarySampleNumber == "{sample_num}" & purpose == "{purpose}"'

                # identify the row with the index sequences in the concatenated fastqFiles + library metadata dataframe
                index_seq_row = metadata_df.query(query_formula)

                # construct a regex style match phrase with the index sequences
                index_1_seq = index_seq_row.index1Sequence.values[0]
                index_2_seq = index_seq_row.index2Sequence.values[0]

                # NOTE: this requires that both index_1 and index_2 be present and correct in the library sheet and the fastq file
                barcode_match = '*' + index_1_seq + '_' + index_2_seq + '*'

                # look for that index in the raw fastq data directory
                for file in os.listdir(fastq_dir_full_path):
                    # if a match is found, enter it in the fastqFile sheet
                    if fnmatch.fnmatch(file, barcode_match):
                        sheet_df.loc[index, 'fastqFileName'] = str(file)

            return sheet_df
                # TODO: fastq_file_name = grep index1seq_index2seq -- if more than one, throw error
                # TODO: enter fastq_file_name into row -- write checks here!


metadata_dir = '/home/chase/code/brentlab/database-files'

# for single sheet input
fastqFile_metadata = '/home/chase/code/brentlab/database-files/fastqFiles/fastq_J.PLAGGENBERG_11.19.19.xlsx'

# for full file input
#fastqFile_metadata = metadata_dir

sequence_dir = '/home/chase/Desktop/sequence'

datadir_keys = ['fastqFiles', 'library']

def main(metadata_dir, datadir_keys, fastqFile_metadata, sequence_dir):
    concat_df = createFastqMetadata(metadata_dir, datadir_keys)
    fastqInput(fastqFile_metadata, concat_df, sequence_dir)


main(metadata_dir, datadir_keys, fastqFile_metadata, sequence_dir)