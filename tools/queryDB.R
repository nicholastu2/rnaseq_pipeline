##################################################################################################################################
# name: queryDB.R
# purpose: create a metadata object from the bio sample datasheets to be inputted to step three of the brent_lab rnaseq pipeline
# input: bio sample database main directory; filepath for script output; json with search terms; name of query
# output: metadata .csv with selected samples
# written by: sanji bhavsar, modified by chase mateusiak
# date included in rnaseq_pipe: 12/19/2019
##################################################################################################################################

# load libraries (TODO: require install.packages if not installed)
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(readxl))
suppressMessages(library(jsonlite))

parseArguments <- function() {
  # parse cmd line arguments
  # Args: flags submitted on cmd line
  # Returns: list of cmd line arguments; error if required flags are not submitted
  
  option_list <- list(
    make_option(c('-i', '--input_directory'), default = NULL,
                help='Path to main directory for database files'),
    make_option(c('-o', '--output'), default = NULL,
                help='Filepath or directory path for newly created metadata sheet'),
    make_option(c('-j', '--json'), default = NULL,
                help='filepath to JSON query'),
    make_option(c('-q', '--query'), default = "query",
                help='some meaningful name for the query')
  )
  
  args <- parse_args(OptionParser(option_list=option_list))
  
  # check if input and output are provided
  if (is.null(args$input_directory) || is.null(args$output) || is.null(args$json)) {
    stop("input, output, and json filepath must be provided.")
  }
  else{
    return(args)
  }
} # end parseArguments()

getFilePaths <- function(bp){
  # create named list of file paths
  # Args: base_path to data directory
  # Returns: named list of all data filepaths in data directory 
  
  fq_files <- list.files(path = bp, pattern = ".astq.*.xlsx", full.names = TRUE, recursive = TRUE)
  lib_files <-  list.files(path = bp, pattern = ".ibrary.*.xlsx", full.names = TRUE, recursive = TRUE)
  s2_files <-  list.files(path = bp, pattern = ".2.DNA.*.xlsx", full.names = TRUE, recursive = TRUE)
  s1_files <-  list.files(path = bp, pattern = ".1.DNA.*.xlsx", full.names = TRUE, recursive = TRUE)
  rna_files <- list.files(path = bp, pattern = ".naSample.*.xlsx", full.names = TRUE, recursive = TRUE)
  bio_files <- list.files(path = bp, pattern = ".ioSample.*.xlsx", full.names = TRUE, recursive = TRUE)
  
  list_names <- c('fq_files', 'lib_files', 's2_files', 's1_files', 'rna_files', 'bio_files')
  
  file_path_list <- list(fq_files, lib_files, s2_files, s1_files, rna_files, bio_files)
  
  names(file_path_list) <- list_names
  
  # halt and print error if no files are found in a given list
  for(i in file_path_list){
    if(length(i) == 0){
      stop(paste("Some files expected to be in the data directory were not found. Check method getFilePath()."))
    }
  }
  
  return(file_path_list)
} # end getFilePaths()

createDB <- function(fpl, output){
  # create joined data frame from data directory
  # Args: fpl is a list of file paths. Must be created from getFilePaths and elements must be in the same order
  # Returns: a complete database of information from data directory, minus entries that do not have fastq paths entered
  
  # create tables from files stored in data directory
  tbl_fastq <- sapply(fpl[[1]], read_excel, simplify=FALSE)  %>%
    bind_rows(.id = "id") %>%
    drop_na(fastqFileName) %>%
    select(-id)
  
  tbl_lib <- sapply(fpl[[2]], read_excel, simplify=FALSE) %>%
    bind_rows(.id = "id") %>% 
    select(-id)
  
  tbl_s2 <- sapply(fpl[[3]], read_excel, simplify=FALSE) %>%
    bind_rows(.id = "id") %>% 
    select(-id)
  
  tbl_s1 <- sapply(fpl[[4]], read_excel, simplify=FALSE) %>%
    bind_rows(.id = "id") %>% 
    select(-id)
  
  tbl_rna <- sapply(fpl[[5]], read_excel, simplify=FALSE) %>%
    bind_rows(.id = "id") %>% 
    select(-id)
  
  tbl_bio <- sapply(fpl[[6]], read_excel, simplify=FALSE) %>%
    bind_rows(.id = "id") %>% 
    select(-id)
  
  # perform all joins
  db <- tbl_fastq %>%
    left_join(tbl_lib) %>%
    left_join(tbl_s2) %>%
    left_join(tbl_s1) %>%
    left_join(tbl_rna) %>%
    left_join(tbl_bio)
  
  # adds vectors holding path to fastq and tsv files
  db <- db %>% mutate(fastqFilePath = paste("sequence/run_", runNumber, "_samples", fastqFileName, sep = ""))
  db <- db %>% mutate(tsvFilePath = gsub(".fastq.gz", "_read_count.tsv", fastqFilePath))
  
  return(db)
} # end createDB

queryDB <- function(df, cols, conds){
  # filter sample dataframe created by createDB for user specified samples
  # Args: a dataframe of all samples in data directory, the column and conditions to filter (supplied in the json via cmd line input)
  # Returns: a filtered dataframe
  
  fp <- map2(cols, conds, function(x, y) quo((!!(as.name(x))) %in% !!y))
  filter(df, !!!fp)
} # end queryDB

writeSummary <- function(query, qname, output){
  
  selectQ <- query %>% select(genotype, strain, inductionDelay, libraryDate, harvestDate, replicate, runNumber, index1Sequence, index2Sequence, fastqFileName, timePoint, floodmedia)
  # print(selectQ)
  summary <- mutate(selectQ,  ST_PIPE = "",
                    ST_TOTAL_READS	= "",
                    ST_ALIGN_PCT	= "",
                    ST_MUT_FOW	= "",
                    ST_RC_FOM	= "",
                    ST_COV_MED	= "",
                    AUTO_AUDIT	= "",
                    MANUAL_AUDIT	= "",
                    USER	= "",
                    NOTE = "")
  
  # write file to csv at specified filepath
  # write_excel_csv(summary, paste(output, qname, "sample_summary.csv", sep = "/"))
  write.csv(x = as.data.frame(summary), file = paste(output, qname, "sample_summary.csv", sep = "/"), row.names = FALSE)
  
  return()
} # end writeSummary

writeLookup <- function(query, qname, output){
  lookupTSV <- tibble(query$tsvFilePath)
  lookupFQ <- tibble(query$fastqFilePath)
  write.table(lookupTSV, paste(output, "/" , qname, "/" , qname, ".expr.lookup.txt", sep = ""), row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)
  write.table(lookupFQ, paste(output, "/" , qname, "/" , qname, ".fastq.lookup.txt", sep = ""), row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)
} # end writeLookup

## Main
# create a joined data.frame from the sample sheets generated in the lab, 
# filter the data.frame for user specified metadata/filepaths, 
# and return a sample_summary sheet 

# store command line arguments as parsed
tryCatch(
parsed <- parseArguments(),
  error = function(a) print('error in parseArguments()')
)

tryCatch(
# creates directory for query specific outputs within output directory
dir.create(file.path(parsed$output, parsed$query, fsep = "/")),
  error = function(b) print('error in main code line ~168. Attempting to create directory in output directory.')
)

# store input directory as base_path
base_path <- as.character(parsed$input_directory)

# get lists of filepaths to relavent data sheets (i.e. fastq files, lib files, etc generated in lab)
file_path_list <- getFilePaths(base_path)

tryCatch(
# create joined database of all relavant sample info and filepaths
db <- createDB(file_path_list, parsed$output),
  error = function(c) print('error in createDB')
)

# create user_query to store user specified filters from json input
json_df <- fromJSON(txt = parsed$json)
user_query <- unnest(json_df, cols = everything())

# filter sample database for user specified items
tryCatch(
filtered_db <- queryDB(db, cols = as.list(names(user_query)), conds = as.list(user_query)),
  error = function(d) print('error in queryDB')
)

# write output files
# write_excel_csv(query, paste(parsed$output, parsed$query, "queriedDB.csv", sep = "/"))
write.csv(x = as.data.frame(filtered_db), file = paste(parsed$output, parsed$query, "queriedDB.csv", sep = "/"), row.names = FALSE)
writeSummary(filtered_db, qname = parsed$query, output = parsed$output)
writeLookup(filtered_db, qname = parsed$query, output = parsed$output)