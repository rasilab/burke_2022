#!/usr/bin/env python
# coding: utf-8

# # Find ribosome density around motifs

# ## Import libraries and set analysis variables

import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
# parameters input to function
parser.add_argument('--read_count_cutoff', nargs='?', default=100, type=int)
parser.add_argument('--read_density_cutoff', nargs='?', default=0.33, type=float)
parser.add_argument('--cds_length_cutoff', nargs='?', default=100, type=int)
parser.add_argument('--l_overhang', nargs='?', default=50, type=int)
parser.add_argument('--r_overhang', nargs='?', default=50, type=int)
parser.add_argument('--transcript_annotation_file', default='/fh/fast/subramaniam_a/db/rasi/genomes/human/hg38/gencode/annotations/gencode.v32.canonical_ccds.parameters.tsv.gz')
parser.add_argument('--motif_annotation_file', default='/fh/fast/subramaniam_a/user/rasi/analysis/bioinformatics/human_codon_usage/tables/ccds_codon_locs.tsv.gz')
parser.add_argument('--ribosome_density_file', default='../data/coverage/GSM2360175.transcripts.bedGraph.gz')
parser.add_argument('--motif_density_file', default='../data/codon_density/GSM2360175.tsv.gz')
parser.add_argument('--average_motif_density_file', default=os.devnull)
args = parser.parse_args()

for k,v in vars(args).items():
    print(k, v)

# class Args(object):
#     pass

# args = Args()

# # parameters input to function
# args.read_count_cutoff = '100'
# args.read_density_cutoff = '0.33' # reads / nt
# args.cds_length_cutoff = '100' # nt
# args.l_overhang = int('50')
# args.r_overhang = int('50')
# args.transcript_annotation_file = '/fh/fast/subramaniam_a/db/rasi/genomes/human/hg38/gencode/annotations/gencode.v32.canonical_ccds.parameters.tsv.gz'
# args.motif_annotation_file = '/fh/fast/subramaniam_a/user/rasi/analysis/bioinformatics/human_codon_usage/tables/20210126_destabilizing_amphipathic_motifs.tsv.gz'
# args.ribosome_density_file = '../data/coverage/GSM1606108.transcripts.bedGraph.gz'
# args.motif_density_file = '../data/vk_type_density//GSM1606108.tsv.gz'
# args.average_motif_density_file = os.devnull


# ## Read in transcript annotations

transcript_annotations = pd.read_table(args.transcript_annotation_file)
transcript_annotations = (
    transcript_annotations
    .set_index('transcript_id')                                       
    .assign(length = lambda x: x.utr5_length + x.cds_length + x.utr3_length)
)
transcript_annotations


# ## Read motifs around which we want to find ribosome density

motif_annotations = pd.read_table(args.motif_annotation_file)
motif_annotations


# ## Read ribosome coverage on transcripts

coverage = pd.read_table(args.ribosome_density_file, names=['transcript_id', 'start', 'end', 'score'])
coverage


# ## Subset to transcripts that have a minimum read count and coverage density

mean_coverage = (
    coverage
    .set_index('transcript_id')
    .join(transcript_annotations)
    .query('start >= utr5_length & end < (utr5_length + cds_length)')
    .reset_index()
    .groupby('transcript_id')[['score']]
    .sum()
    .join(transcript_annotations)
    .assign(cds_read_density = lambda x: x.score / x.cds_length)
    .query(f'score > {args.read_count_cutoff}')
    .query(f'cds_read_density > {args.read_density_cutoff}')
    .query(f'cds_length > {args.cds_length_cutoff}')
    .sort_values('cds_read_density', ascending=False)
    .reset_index()
    .rename({'score': 'cds_read_count'}, axis=1)
    .loc[:, ['transcript_id', 'cds_read_density', 'cds_read_count']]
)

mean_coverage


# ## Subset coverage to transcripts that pass cutoff

subset_coverage = (
    coverage
    .query('score > 0')
    .merge(mean_coverage, on = 'transcript_id', how='inner')
)

subset_coverage


# # Convert coverage to numpy array

coverage_array = dict()
for row in mean_coverage.iterrows():
    coverage_array[row[1]['transcript_id']] = np.zeros(
        transcript_annotations.loc[row[1]['transcript_id'], 'length'].item() + args.l_overhang + args.r_overhang)
    coverage_array[row[1]['transcript_id']][:args.l_overhang] = np.nan
    coverage_array[row[1]['transcript_id']][-args.r_overhang:] = np.nan
    
for row in subset_coverage.iterrows():
    tx = row[1]['transcript_id']
    start =  row[1]['start'] + args.l_overhang
    end =  row[1]['end'] + args.r_overhang
    score =  row[1]['score']
    coverage_array[tx][start:end] = score


# ## Find motifs only in transcripts that pass cutoff and expand motif on either side by `l_overhang` and `r_overhang`

subset_motif_annotations = (
    motif_annotations
    .merge(mean_coverage, on = 'transcript_id', how='inner')
    .assign(start = lambda x: x['loc'])
    .assign(end = lambda x: x['loc'] + args.l_overhang + args.r_overhang)
)

subset_motif_annotations


# ## Create empty array to hold motif coverage

coverage_col_names = list(range(-args.l_overhang,0)) + list(range(1, args.r_overhang + 1))
motif_coverage = np.empty((len(subset_motif_annotations),args.l_overhang + args.r_overhang))
motif_coverage.fill(None)
motif_coverage


# ## Calculate coverage around each motif

motif_n = 0
for row in subset_motif_annotations.iterrows():
    tx = row[1]['transcript_id']
    start =  row[1]['start']
    end =  row[1]['end'] 
    motif_coverage[motif_n, :] = coverage_array[tx][start:end]
    motif_n += 1
    if motif_n % 1e5 == 0:
        print(motif_n)


# Convert motif coverage to a DataFrame for joining with motif annotations

motif_coverage_df = pd.DataFrame(motif_coverage, dtype=pd.Int64Dtype)
motif_coverage_df.columns = coverage_col_names
motif_coverage_df = pd.concat([subset_motif_annotations, motif_coverage_df], axis=1)
motif_coverage_df


# ## Write unnormalized coverage around each motif to output

(
    motif_coverage_df
    .assign(cds_read_density = lambda x: np.round(x.cds_read_density, 3))
    .sort_values('cds_read_density', ascending=False)
    .to_csv(args.motif_density_file, sep='\t', index=False)
)


# ## Calculate average coverage across all occurrences of motif after normalizing by CDS read density

average_coverage = {k:np.nan for k in subset_motif_annotations.motif.unique()}

for motif in average_coverage:
    idx = subset_motif_annotations['motif'] == motif
    a = motif_coverage[idx,:] 
    b = subset_motif_annotations.loc[idx, 'cds_read_density'][:, None]
    average_coverage[motif] = np.nanmean( a / b, axis=0)


average_coverage = (
    pd.DataFrame(average_coverage)
    .transpose()
    .rename({n:coverage_col_names[n] for n in range(100)}, axis=1)
    .apply(lambda x: np.round(x, 4))
)
average_coverage


# ## Write average motif density to output file

average_coverage.to_csv(args.average_motif_density_file, "\t", index_label='motif')

