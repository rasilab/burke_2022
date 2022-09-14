#' ---
#' title: Calculate counts from transcript coverage
#' author: Arvind Rasi Subramaniam
#' date: 27 Jan 2021
#' ---
#' 
#' **Edit this Rscript only in the accompanying .Rmd file with same name and
#' export by running the last cell in the .Rmd file.**
#'   
#' ## Import libraries and arguments
## --------------------------------------------------------------------------------------------------------------------------------------------------
# for UCSC tracks
library(rtracklayer)
# for string concatenation
library(glue)
# for tab data work
library(tidyverse)
# for tidy operations on Bioconductor objects
library(plyranges)

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_file <- args[2]

#' 
#' ## Read coverage and write summed coverage to output
## --------------------------------------------------------------------------------------------------------------------------------------------------
coverage <- import.bedGraph(input_file) %>% 
  as_tibble() %>% 
  group_by(seqnames) %>% 
  summarize(score = sum(score)) %>% 
  ungroup() %>% 
  rename(transcript_id = seqnames) %>% 
  arrange(-score) %>% 
  write_tsv(output_file) %>% 
  print()

#' 
