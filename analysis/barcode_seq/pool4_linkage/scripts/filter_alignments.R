#' ---
#' title: Find barcodes for each insert
#' author: Arvind Rasi Subramaniam
#' date: 26 Mar 2021
#' ---
#' 
#' **Edit this Rscript only in the accompanying .Rmd file with same name and
#' export by running the last cell in the .Rmd file.**
#' 

#' 
#' # Load libraries and define analysis-specific parameters
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
input_file <- args[1]
output_file <- args[2]
input_file <- "../data/alignments/stall4control4_linkage.bam"
output_file <- "../data/filtered_alignments/stall4control4_linkage.tsv.gz"
print(input_file)

library(Biostrings)
library(GenomicAlignments)
library(plyranges)
library(tidyverse)

#' 
#' # Read alignments
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
bamfile <- BamFile(input_file)
# extract the number of mismatches and total edits
param <- ScanBamParam(
  # what = scanBamWhat(), 
  what = c('qname', 'flag'),
  # no secondary alignments
  flag = scanBamFlag(isSecondaryAlignment = F),
  # extract number of mismatches
  tag = c("XM"), 
  # # limit number of mismatches to 1, 
  # tagFilter = list("XM" = c(1)),
  # include only snps; exclude indels
  simpleCigar = T,
  # skip reads that have >1% chance of mapping to another ref
  mapqFilter = 20
)
alns <- readGAlignments(bamfile, param = param) 

#' 
#' 
#' # Write read-ref pairs to output
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
alns %>% 
  GRanges() %>% 
  as_tibble() %>% 
  select(seqnames, qname, width) %>%
  rename(ref = seqnames) %>% 
  write_tsv(output_file)
  # group_by(XM) %>% 
  # summarize(n_mismatch = dplyr::n()) %>% 
  # print()

#' 
#' 
