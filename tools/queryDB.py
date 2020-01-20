#!/usr/bin/env python

# name: queryDB
# purpose: parse rnaseq metadata
# input: bio sample database main directory; filepath for script output; json with search terms; name of query
# output: four .csv : expr.lookup.txt has paths to the count matricies;
#                     fastq.lookup.txt has paths to the fastqs;
#                     queriedDB.csv is the complete database joined from input directory;
#                     sample_summary.csv is the filtered table for input into rnaseq_pipe
# written by: implemented in python by chase mateusiak(chase.mateusiak@gmail.com), based on sanji bhavsar code in R.
#             See https://github.com/sanjibhavsar/rnaseq_database for original work
# date included in rnaseq_pipe: 1/20/2020
# See the end of the script for a description of the environment used to write/test

import os
import pandas as pd
import re
import sys

# list of subdirectories of datadir to be used as keys of dictionary. When these files are searched, only .csv and .xlsx
# are listed for concatenating. The search through these directories is not recursive.
datadir_keys = ['fastqFiles', 'library', 's2cDNASample', 's1cDNASample', 'rnaSample', 'bioSample']

datadir = '/home/chase/code/brentlab/database-files'
query = pd.read_json('/home/chase/Desktop/json_query3.json', typ='series')

def main(datadir, query):
    datadir_dict = getFilePaths(datadir)

    combined_df = createDB(datadir_dict)

    queried_df = queryDB(combined_df, query)

    queried_df.to_csv("~/Desktop/queryDB.csv")

def parseArgs(argv):
    parser = argeparse.ArgumentParser()
    parser.add_argument('-d', '--database', required = True,
                        help = 'topmost directory of metadata database. On cluster, /scratch/mblab/database-files. \
                                If using the rnaseq_pipeline module, you may use $METADATA. Do make sure that $METADATA \
                                is up to date by running git pull in the directory')
    parser.add_argument('-j' '--json', required = True,
                        help = 'path to json file used to parse metadata. See ')
    parser.add_argument('-o,' '--output', required = True,
                        help = 'filepath to directory to intended queryDB output')
    parser.add_argument('-q', '--query', required = True,
                        help = 'name of the query. This will be used to name files deposited in the specified output \
                        directory')

def getFilePaths(datadir, datadir_keys = datadir_keys):
    # create dictionary of filepaths to the various types of metadata sheets
    # Args: base_path to data directory
    # Returns: dictionary {metadata-database/subdir: [filepaths]} i.e. {bioSample: [filepaths], experimentDesign:[filepaths]}

    datadir_dict = {}

    # associate each key (relevant subdirs of datadir) with a list of of files (complete path) in each key directory
    for key in datadir_keys:
        dir_path = os.path.join(datadir, key)
        subdir_files = []
        for filename in os.listdir(os.path.join(datadir,key)):
            if(not filename.startswith('._') and re.search('(xlsx|csv)',filename)): # no hidden files or sub directories
                subdir_files.append(os.path.join(dir_path,filename))

        # test whether any of the key subdirectories of datadir are empty, throw error if so
        try:
            1/len(subdir_files)
        except ZeroDivisionError:
            print("No files found in %s. These files are necessary to creating the sample_summary." % os.path.join(datadir,key))
            exit(0)
        else:
            pass
        datadir_dict.setdefault(key, []).append(subdir_files)

    return datadir_dict

def checkCSV(file):
    # test whether a given file is a .csv or .xlsx
    if re.search('\.csv', file):
        return True
    else:
        return False

def concatMetadata(metadata_sheet_list):
    # creates concatenatd dataframe from all files in a given list of paths (files of a certain subdirectory of user inputted data directory
    # Args: a list of filepaths all pointing to sheets of the same structure
    # Returns: all of the sheets concatenated vertically

    # concatenate (row_bind) all tables in a given directory's list of filepaths
    if checkCSV(metadata_sheet_list[0][0]):
        concatenated_df = pd.read_csv(metadata_sheet_list[0][0])
    else:
        concatenated_df = pd.read_excel(metadata_sheet_list[0][0])
    for path in metadata_sheet_list[0][1:]:
        if checkCSV(path):
            next_df = pd.read_csv(path)
            concatenated_df = concatenated_df.append(next_df)
        else:
            next_df = pd.read_excel(path)
            concatenated_df = concatenated_df.append(next_df)
    return concatenated_df

def createDB(datadir_dict, datadir_keys = datadir_keys, drop_fastq_na = True):
    # create joined data frame from data directory
    # Args: fpl is a list of file paths. Must be created from getFilePaths and elements must be in the same order
    # Returns: a complete database of information from data directory, minus entries that do not have fastq paths entered

    concat_dict = {}

    # return dictionary of structure {datadir/subdir: concatenated_table_of_all_files_in_datadir/subdir}
    for key in datadir_dict:
        concat_dict[key] = concatMetadata(datadir_dict[key])

    # drop rows with na in column fastqFileName of table fastqFiles
    if drop_fastq_na:
        concat_dict['fastqFiles'].dropna(subset=['fastqFileName'], inplace=True)

    # create list of shared columns between successive pairs of keys in datadir_keys (global variable at top). i.e. the
    # columns which are shared between the sheets in fastqFiles and Library
    key_cols = []
    for i in range(len(datadir_keys)-1):
        key_cols.append(concat_dict[datadir_keys[i]].keys().intersection(concat_dict[datadir_keys[i+1]].keys()))

    # merge the first two sets of data, the concatenated fastqFiles and Library sheets
    merged_df = pd.merge(concat_dict[datadir_keys[0]], concat_dict[datadir_keys[1]], how='left', on=list(key_cols[0]))

    # merged_df.to_csv('~/Desktop/merged_df_1.csv') # debugging tool. this and the print line below will print out each intermediate merged_df
    # merge the subsequent sheets on the columns identified in key_cols

    for i in range(1,len(datadir_keys)-1):
        merged_df = pd.merge(merged_df, concat_dict[datadir_keys[i+1]], how='left', on=list(key_cols[i]))

        # name = "~/Desktop/merged_df_"+str(i+1)+".csv" # debugging tool
        # merged_df.to_csv(name)

    merged_df = merged_df.astype('str')
    return merged_df

def queryDB(df, query):
    # filters combined_df on user-inputted columns/values
    # Args: query (user input json at cmd line)
    # Return: a filtered df

    #TODO: make query case insensitive

    # begin a string to store the query formula
    fltr_str='('
    # loop through columns in json query (i.e. 'timePoint' and 'treatment')S
    for key in query.index:
        # loop through filters for each column in json query (i.e. [15, -1] and [EtoH])
        for value in query[key]:
            # for each value, create a formula in the form '(column == "value") |
            fltr_str = fltr_str + '({} == "{}") | '.format(str(key), str(value))
        # once finished looping through values of given column, close the paren opened in the instantiation of fltr_str
        # and add a & ( to capture the next column on which to filter
        fltr_str = fltr_str[:-3] + ') & ('
    # once out of both loops, split the string on the final &, retain the substring before the split, and eliminate whitespace
    fltr_str = fltr_str.rsplit('&',1)[0].strip()
    # use the fltr_str formula to filter the dataframe
    df = df.query(fltr_str)
    return df

if __name__ == '__main__':
	main(sys.argv)

# environment used to write/test
# platform: linux-64
# @EXPLICIT
# https://conda.anaconda.org/conda-forge/linux-64/_libgcc_mutex-0.1-conda_forge.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/ca-certificates-2019.11.28-hecc5488_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/ld_impl_linux-64-2.33.1-h53a641e_7.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/libgfortran-ng-7.3.0-hdf63c60_4.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/libstdcxx-ng-9.2.0-hdf63c60_2.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/libgomp-9.2.0-h24d8f2e_2.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/_openmp_mutex-4.5-0_gnu.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/libgcc-ng-9.2.0-h24d8f2e_2.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/libffi-3.2.1-he1b5a44_1006.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/libopenblas-0.3.7-h5ec1e0e_6.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/ncurses-6.1-hf484d3e_1002.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/openssl-1.1.1d-h516909a_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/xz-5.2.4-h14c3975_1001.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/yaml-0.2.2-h516909a_1.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/zlib-1.2.11-h516909a_1006.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/libblas-3.8.0-14_openblas.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/readline-8.0-hf8c457e_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/tk-8.6.10-hed695b0_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/libcblas-3.8.0-14_openblas.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/liblapack-3.8.0-14_openblas.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/sqlite-3.30.1-hcee41ef_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/python-3.7.6-h357f687_2.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/certifi-2019.11.28-py37_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/chardet-3.0.4-py37_1003.tar.bz2
# https://conda.anaconda.org/conda-forge/noarch/et_xmlfile-1.0.1-py_1001.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/idna-2.8-py37_1000.tar.bz2
# https://conda.anaconda.org/conda-forge/noarch/jdcal-1.4.1-py_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/numpy-1.17.3-py37h95a1406_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/pycosat-0.6.3-py37h516909a_1002.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/pycparser-2.19-py37_1.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/pysocks-1.7.1-py37_0.tar.bz2
# https://conda.anaconda.org/conda-forge/noarch/pytz-2019.3-py_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/ruamel_yaml-0.15.80-py37h516909a_1000.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/six-1.13.0-py37_0.tar.bz2
# https://conda.anaconda.org/conda-forge/noarch/tqdm-4.41.1-py_0.tar.bz2
# https://conda.anaconda.org/conda-forge/noarch/xlrd-1.2.0-py_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/cffi-1.13.2-py37h8022711_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/conda-package-handling-1.6.0-py37h516909a_1.tar.bz2
# https://conda.anaconda.org/conda-forge/noarch/openpyxl-3.0.1-py_0.tar.bz2
# https://conda.anaconda.org/conda-forge/noarch/python-dateutil-2.8.1-py_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/setuptools-44.0.0-py37_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/cryptography-2.8-py37h72c5cf5_1.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/pandas-0.25.3-py37hb3f55d8_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/wheel-0.33.6-py37_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/pip-19.3.1-py37_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/pyopenssl-19.1.0-py37_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/urllib3-1.25.7-py37_0.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/requests-2.22.0-py37_1.tar.bz2
# https://conda.anaconda.org/conda-forge/linux-64/conda-4.8.0-py37_1.tar.bz2
