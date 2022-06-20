#' @title Create bowtie reference from oligo pool sequences 
#' 
#' @author  Arvind R. Subramaniam
#' Date: 30 Oct 2020
#' 
#'
#' @param oligo_pool_file Name of file in TSV or TXT format
#' @param output_file Name of file in FASTA format

library(Biostrings)
library(tidyverse)



args <- commandArgs(trailingOnly = T)
oligo_pool_file <- args[1]
output_file <- args[2]

readLines(oligo_pool_file) %>% 
  enframe("insert_num", "oligo_seq") %>% 
  mutate(insert_num = 0:(n()-1)) %>% 
  select(insert_num, oligo_seq) %>%
  mutate(oligo_seq = toupper(oligo_seq)) %>% 
  mutate(insert = str_extract(oligo_seq, 
                                 "(?<=CTAAGGCCATCACCAAGTACACCAGCGCT)[ATCG]+(?=AGCAAGGGCGAGGAGCTGTTCACC)")) %>%
  select(insert_num, oligo_seq, insert) %>%
  write_tsv("../annotations/insert_annotations.tsv.gz") %>%
  mutate(oligo_seq = str_extract(oligo_seq, 
                                 "CTAAGGCCATCACCAAGTACACCAGCGCT[ATCG]+AGCAAGGGCGAGGAGCTGTTCACC")) %>%
  select(insert_num, oligo_seq) %>%
  deframe() %>%
  DNAStringSet() %>%
  reverseComplement() %>% 
  writeXStringSet(output_file)
