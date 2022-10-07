#' ---
#' title: Filter barcodes to remove ones aligning to multiple inserts or second barcode
#' author: Arvind Rasi Subramaniam
#' date: Mar 26, 2021
#' ---
#' 
#' **Edit this Rscript only in the accompanying .Rmd file with same name and
#' export by running the last cell in the .Rmd file.**
#' 

#' 
#' # Load libraries and define analysis-specific parameters
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
barcode1_alignment_file <- args[1]
barcode_insert_file <- args[2]
output_file <- args[3]
# barcode1_alignment_file <- '../data/ref_vs_ref_alignments/endo12k/alignment_barcode1.bam'
# barcode_insert_file <- '../data/insert_barcode_counts/endo12k.tsv.gz'
# output_file <- '../data/filtered_barcodes/endo12k.tsv.gz'

library(Biostrings)
library(GenomicAlignments)
library(plyranges)
library(tidyverse)

#' 
#' # Read barcodes that we spike-into all samples
#' 
#' 
#' # Read insert-barcode pair counts 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
insert_barcodes <- read_tsv(barcode_insert_file) %>% 
  print()

#' 
#' 
#' # How many barcode_1 have multiple inserts?
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
many_to_one_barcode_combinations <- insert_barcodes %>% 
  group_by(barcode_1) %>% 
  mutate(n1 = dplyr::n()) %>% 
  ungroup() %>% 
  filter(n1 > 1) %>% 
  print()

#' 
#' 
#' # Fields to read from BAM file
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# extract the number of mismatches and total edits
param <- ScanBamParam(
  # what = scanBamWhat(),
  what = c("qname", "flag"),
  # extract number of mismatches
  tag = c("XM"), 
  # include only snps; exclude indels
  simpleCigar = T
)

#' 
#' 
#' # Read barcode vs barcode alignments for barcodes 1
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
bamfile1 <- BamFile(barcode1_alignment_file)
alns1 <- readGAlignments(bamfile1, param = param) %>% 
  as_tibble() %>% 
  mutate(rname = as.character(seqnames)) %>% 
  select(rname, qname, flag, XM) %>% 
  type_convert() %>% 
  print()

#' 
#' 
#' # Find barcode_1 that are linked to distinct insert or might be sequencing errors
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
exclude1 <- alns1 %>% 
  filter(rname != qname) %>%
  left_join(select(insert_barcodes, insert_num, barcode_num, read_count), by = c("rname" = "barcode_num")) %>%
  rename(rinsert = insert_num, rcount = read_count) %>%
  right_join(select(insert_barcodes, insert_num, barcode_num, read_count), by = c("qname" = "barcode_num")) %>%
  rename(qinsert = insert_num, qcount = read_count) %>%
  # this exludes:
  # 1. barcodes that map to two distinct inserts
  # 2. barcodes that got lower count than another homologous barcode with same insert
  filter(!(qinsert == rinsert & qcount > rcount)) %>%
  arrange(qname) %>% 
  distinct(qname) %>%
  print()

#' 
#' # Write barcodes that do not clash to output
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
filtered_barcodes <- insert_barcodes %>% 
  anti_join(select(exclude1, qname), by = c("barcode_num" = "qname")) %>%
  anti_join(select(many_to_one_barcode_combinations, barcode_num), by = "barcode_num") %>%
  select(insert_num, barcode_num, barcode_1, read_count) %>%
  arrange(desc(read_count)) %>%
  mutate(barcode_num = 1:dplyr::n()) %>%
  write_tsv(output_file) %>%
  print()

#' 
#' 
