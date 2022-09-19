#' ---
#' title: Create bowtie reference
#' author: Arvind Rasi Subramaniam
#' date: 27 Mar 2021
#' ---
#' 
#' **Edit this Rscript only in the accompanying .Rmd file with same name and
#' export by running the last cell in the .Rmd file.**
#' 

#' 
#' # Load libraries and define analysis-specific parameters
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
linkage_ref_file <- args[1]
barcode1_fasta_file <- args[2]
# linkage_ref_file <- "../annotations/insert_barcode_linkages/didi8k.tsv.gz"
# barcode1_fasta_file <- "../data/bowtie_references/didi8k/reference_barcode_1.fasta"
# barcode2_fasta_file <- "../data/bowtie_references/didi8k/reference_barcode_2.fasta"

library(Biostrings)
library(tidyverse)
library(GenomicRanges)

#' 
#' # Create reference file for barcode 1
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
linkage_ref_file %>% 
  read_tsv() %>% 
  select(barcode_num, barcode_1) %>%
  deframe() %>%
  DNAStringSet() %>%
  writeXStringSet(barcode1_fasta_file) %>% 
  print()

#' 
