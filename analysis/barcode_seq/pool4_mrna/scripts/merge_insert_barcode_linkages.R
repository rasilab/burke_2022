#' ---
#' title: Merge insert barcode linkages from different plasmid pools
#' author: Arvind Rasi Subramaniam
#' date: 3 Sep 2021
#' ---
#' 
#' **Edit this Rscript only in the accompanying .Rmd file with same name and
#' export by running the last cell in the .Rmd file.**
#' 

#' 
#' # Load libraries and define analysis-specific parameters
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

#' 
#' # Read insert barcode linkages and merge them
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
linkage <- list.files("../../pool4_linkage/data/filtered_barcodes/", pattern = "stall4control4_", full.names = T) %>% 
  enframe("sno", "file") %>% 
  mutate(data = map(file, read_tsv)) %>% 
  mutate(pool = str_extract(file, "[:digit:](?=.tsv.gz)")) %>% 
  select(-sno, -file) %>% 
  unnest(data) %>% 
  mutate(barcode_num = seq(1,dplyr::n())) %>% 
  write_tsv("../annotations/insert_barcode_linkages/stall4control4.tsv.gz") %>% 
  print()

#' 
#' 
#' 
