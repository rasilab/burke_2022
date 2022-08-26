#' ---
#' title: Create references for aligning 12K and 8K libraries
#' author: Arvind Rasi Subramaniam
#' date: 22 Dec 2020
#' ---
#' 
#' **Edit this Rscript only in the accompanying .Rmd file with same name and
#' export by running the last cell in the .Rmd file.**
#' 

#' 
#' # Load libraries and define analysis-specific parameters
## --------------------------------------------------------------------------------------------------------------------------------------------------
library_annotations_files <- c(human = "/fh/fast/subramaniam_a/user/rasi/analysis/cloningdesign/20210216_endogenous_vk_oligo_pool_design/tables/20210216_vk_type_endogenous_motifs_for_cloning.tsv",
                               viral = "/fh/fast/subramaniam_a/user/rasi/analysis/cloningdesign/20210216_endogenous_vk_oligo_pool_design/tables/20210217_viral_rg4_motifs_for_cloning.tsv")
output_file <- "../data/bowtie2_reference/endo12k.fasta"

# from pHPHS286
flanking_5 <- "GGTGAACAGCTCCTCGCCCTTGCT"
flanking_3 <- "agcgctggtgtacttggtgatgGCCTTAG"

# maximum insert size
max_insert_size = 48 

# how many 5' nt to include in ref
adapter_5_size = 0

# how many 3' nt to include in ref
adapter_3_size = 0

# trim insert + flanking region starting at this location
seq_start <- nchar(flanking_5) - (adapter_5_size - 1)

# trim insert  + flanking region at this length
seq_length <- adapter_5_size + max_insert_size + adapter_3_size

library(Biostrings)
library(tidyverse)
library(GenomicRanges)

#' 
#' # Use values below for trimming the reads during bowtie2 alignment
#' 
## --------------------------------------------------------------------------------------------------------------------------------------------------
seq_start
seq_length

#' 
#' # Create reference sequences for alignment
## --------------------------------------------------------------------------------------------------------------------------------------------------
library_annotations_files %>% 
  enframe("group", "file") %>% 
  mutate(data = map(file, read_tsv)) %>%
  select(-file) %>%
  unnest() %>%
  arrange(desc(group)) %>%
  mutate(insert_num = 0:(dplyr::n()-1)) %>%
  arrange(-insert_num) %>% 
  write_tsv("../annotations/insert_annotations/endo12k.tsv") %>%
  arrange(insert_num) %>% 
  mutate(seq = subseq(seq, end = width(seq), width = seq_length)) %>%
  rename(insert_seq = seq) %>%
  select(insert_num, insert_seq) %>%
  deframe() %>%
  DNAStringSet() %>%
  writeXStringSet(output_file)
  # print()

#' 
#' 
