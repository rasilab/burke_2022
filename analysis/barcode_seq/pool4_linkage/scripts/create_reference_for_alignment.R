#' ---
#' title: Create references for aligning 12K and 8K libraries
#' author: Arvind Rasi Subramaniam
#' date: 22 Dec 2020
#' ---

#' 
#' # Load libraries and define analysis-specific parameters
## --------------------------------------------------------------------------------------------------------------------------------------------------
library(Biostrings)
library(tidyverse)
library(GenomicRanges)

library_annotations_files <- c("../annotations/insert_annotations.tsv")
output_file <- "../data/bowtie2_reference/stall4control4_linkage.fasta"

#' 
#' # Create reference sequences for alignment
## --------------------------------------------------------------------------------------------------------------------------------------------------
library_annotations_files %>% 
  enframe("group", "file") %>% 
  mutate(data = map(file, read_tsv)) %>%
  select(-file) %>%
  unnest() %>%
  write_tsv("../annotations/insert_annotations/stall4control4_linkage.tsv") %>%
  arrange(insert_num) %>% 
  select(insert_num, insert_seq) %>%
  deframe() %>%
  DNAStringSet() %>%
  writeXStringSet(output_file)