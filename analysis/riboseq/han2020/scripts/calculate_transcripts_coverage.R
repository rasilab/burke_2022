#' ---
#' title: Calculate coverage from transcript alignments
#' author: Arvind Rasi Subramaniam
#' date: 29 Oct 2020
#' ---
#' 
#' **Edit this Rscript only in the accompanying .Rmd file with same name and
#' export by running the last cell in the .Rmd file.**
#' 
#' ## Load libraries
## -------------------------------------------------------------------------------------------------------------------------------------------------
# for UCSC tracks
library(rtracklayer)
# for handling htseq alignments
library(GenomicAlignments)
# for string concatenation
library(glue)
# for tab data work
library(tidyverse)
# for tidy operations on Bioconductor objects
library(plyranges)

#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
output_file <- args[1]
trim_file <- args[2]
study_name <- args[3]
input_files <- args[4:length(args)]

# output_file <- '../data/coverage/GSM2360179.transcripts.bedGraph.gz'
# trim_file <- '../annotations/p_site_trim_distances.csv'
# study_name <- 'khajuria2018'
# input_files <- c('../data/alignments/SRR4450327.transcripts.bam')

#' 
#' 
#' ## Read trim parameters file containing allowed read lengths and trimming distance to P-site
## -------------------------------------------------------------------------------------------------------------------------------------------------
trim_lengths <- read_csv(trim_file) %>% 
  filter(study == study_name)

left_trim <- setNames(trim_lengths$left_trim, as.character(trim_lengths$read_length))
right_trim <- setNames(trim_lengths$right_trim, as.character(trim_lengths$read_length))
allowed_lengths <- trim_lengths$read_length
remaining_width <- unlist(map2(right_trim, left_trim, function(x, y) x - y + 1))

#' 
#' ## Read in alignments
## -------------------------------------------------------------------------------------------------------------------------------------------------
# retrieve alignments 
message(glue("Reading alignments."))
# skip secondary alignments
param <- ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE), what=c("flag"))

alns <- sapply(input_files, . %>% BamFile %>% readGAlignments(param=param))
alns <- unlist(GAlignmentsList(alns), use.names = F)
print(alns)

#' 
#' ## Trim alignments
## -------------------------------------------------------------------------------------------------------------------------------------------------
# length of query
alnLength <- qwidth(alns)
# how much to trim on left side
alnLeftTrim <- left_trim[as.character(alnLength)]
# Width of remaining alignment
alnWidth <- remaining_width[as.character(alnLength)]
# Subset to alignments that have positive trimmed length
alnGood <- alnLength %in% allowed_lengths
if (length(alns[alnGood]) > 0) {
  alnStart <- ifelse(
    strand(alns[alnGood]) == "+",
    # trim by left.trim if aln is on + strand
    alnLeftTrim[alnGood] + 1,
    # trim by aln_length - left.trim if aln is on - strand
    alnLength[alnGood] - alnLeftTrim[alnGood])
  psites <- qnarrow(alns[alnGood],
                    start = alnStart,
                    width = alnWidth[alnGood])
} else {
  # return empty range if no aln pass the above checks
  psites <- GRanges()
}

#' 
#' ## Calculate coverage
## -------------------------------------------------------------------------------------------------------------------------------------------------
cvg <- GRanges(coverage(psites)) %>% 
  arrange(seqnames, start) %>%
  print()

#'   
#' ## Export coverage
## -------------------------------------------------------------------------------------------------------------------------------------------------
export(cvg,
       con=output_file,
       format="bedGraph")
message(glue("Exported {output_file}.\n"))

#' 
