# Extract unique CCDS for each gene from Gencode annotations ----

## Load libraries ----
# %% 
library(rtracklayer)
library(Biostrings)
library(plyranges)
library(tidyverse)
# %%

## Read genome ----
# %%
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
# %%

## Read Gencode annotations ----
# %%
gencode_annotations <- import.gff3("../data/gencode/gencode.v32.annotation.gff3.gz")
# %%

#'
#' ## Get one transcript per CCDS
# %%
## -----------------------------------------------------------------------------
ccds <- gencode_annotations %>%
  filter(type == "CDS") %>%
  as_tibble() %>%
  # select only CCDS
  filter(!is.na(ccdsid) | (transcript_support_level == 1)) %>%
  # skip transcripts in chrY if they are already in X
  group_by(transcript_id) %>%
  mutate(n_chr = n_distinct(seqnames)) %>%
  ungroup() %>%
  filter(!(n_chr > 1 & seqnames == "chrY")) %>%
  select(-n_chr) %>%
  mutate(length = abs(start - end)) %>%
  group_by(gene_id, ccdsid) %>%
  arrange(transcript_id) %>%
  # select lowest transcript_id for each CCDS
  filter(transcript_id == min(transcript_id)) %>%
  mutate(ccds_length = sum(length)) %>%
  ungroup() %>%
  group_by(gene_id) %>%
  # select longest CCDS for each gene
  filter(ccds_length == max(ccds_length)) %>%
  # pick only the first CCDS after this
  filter(ccdsid == min(ccdsid) | is.na(ccdsid)) %>%
  ungroup() %>%
  GRanges() %>%
  print()
#%%

#' ## Get CCDS annotations
## -----------------------------------------------------------------------------
ccds_annotations <- ccds %>%
  mcols() %>%
  as_tibble() %>%
  select(transcript_id, gene_id, gene_name, ccdsid) %>%
  group_by(transcript_id) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  arrange(transcript_id) %>%
  print()

#'
#' ## Get nucleotide sequence of all CCDS
## -----------------------------------------------------------------------------
ccds_seqs <- ccds %>%
  arrange(transcript_id) %>%
  split(.$transcript_id) %>%
  GenomicFeatures::extractTranscriptSeqs(genome, .) %>%
  print()

# write to file
writeXStringSet(ccds_seqs, file = "../data/gencode/gencode.v32.canonical_ccds.fa", format = "fasta")

#'
#'
#' ## Get annotations of all CCDS transcripts (including 5'UTR and 3'UTR)
## -----------------------------------------------------------------------------
ccds_tx <- gencode_annotations %>%
  as_tibble() %>%
  filter(type %in% c("CDS", "five_prime_UTR", "three_prime_UTR")) %>%
  filter(transcript_id %in% ccds_annotations$transcript_id) %>%
  # skip transcripts in chrY if they are already in X
  group_by(transcript_id) %>%
  mutate(n_chr = n_distinct(seqnames)) %>%
  ungroup() %>%
  filter(!(n_chr > 1 & seqnames == "chrY")) %>%
  select(-n_chr) %>%
  # arrange exons 5' to 3'
  group_by(transcript_id) %>%
  mutate(order = if_else(strand == "+", start, -start)) %>%
  arrange(order) %>%
  ungroup() %>%
  select(-order) %>%
  write_tsv("../data/gencode/gencode.v32.canonical_ccds_transcripts.annotations.20210126.tsv.gz") %>%
  GRanges()

#'
#'
#' ## Calculate CDS, UTR lengths and write to file
## -----------------------------------------------------------------------------
ccds_tx %>%
  as_tibble() %>%
  group_by(transcript_id) %>%
  mutate(
    utr5_length = if_else(type == "five_prime_UTR", width, as.integer(0)),
    cds_length = if_else(type == "CDS", width, as.integer(0)),
    utr3_length = if_else(type == "three_prime_UTR", width, as.integer(0))
  ) %>%
  summarize(utr5_length = sum(utr5_length), cds_length = sum(cds_length), utr3_length = sum(utr3_length)) %>%
  inner_join(ccds_annotations, by = "transcript_id") %>%
  write_tsv("../data/gencode/gencode.v32.canonical_ccds.parameters.tsv.gz")

#'
#' ## Get nucleotide sequence of all CCDS transcripts and write to file
## -----------------------------------------------------------------------------
ccds_tx_seqs <- ccds_tx %>%
  split(.$transcript_id) %>%
  GenomicFeatures::extractTranscriptSeqs(genome, .) %>%
  print()

writeXStringSet(ccds_tx_seqs, file = "../data/gencode/gencode.v32.canonical_ccds_tx.fa", format = "fasta")

#'