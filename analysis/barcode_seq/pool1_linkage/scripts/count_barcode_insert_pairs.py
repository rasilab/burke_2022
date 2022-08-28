#!/usr/bin/env python
# coding: utf-8

# # Read filtered alignments and find and count corresponding barcodes
# 
# We include only insert-barcode pairs that are supported 10 or more reads.

import pandas as pd
import numpy as np
import sys
import plotnine as pn
from Bio.SeqIO.QualityIO import FastqGeneralIterator


# ## Define analysis-specific parameters

r1_fastq_file = sys.argv[1]
alignment_file = sys.argv[2]
output_file = sys.argv[3]
barcode1_read = int(sys.argv[4])
barcode1_start = int(sys.argv[5])
barcode1_length = sys.argv[6]
# barcode2_read = int(sys.argv[7])
# barcode2_start = int(sys.argv[8])
# barcode2_length = sys.argv[9]

# r1_fastq_file = '../data/fastq/139P7_S3_R1_001.fastq.gz'
# alignment_file = '../data/filtered_alignments/endo12k.tsv.gz'
# output_file = '../data/insert_barcode_counts/endo12k.tsv.gz'
# barcode1_read = 1
# barcode1_start = 0 
# barcode1_length = '24'
# barcode2_read = 1
# barcode2_start = 67
# barcode2_length = '19'


barcode1_end = barcode1_start + int(barcode1_length)
# barcode2_end = barcode2_start + int(barcode2_length)


# ## Read in the filtered alignments table

aln = pd.read_table(alignment_file).sort_values(by = 'qname').set_index('qname')


# ## Parse dual barcodes from Read 1 or Read 2 for each filtered alignment
# 
# Note that the params `barcode(1|2)_(read|start|length)` specify which read and what location the barcodes are present.

f1 = FastqGeneralIterator(open(r1_fastq_file, 'rt'))
# r2_fastq_file = r1_fastq_file.replace('R1', 'R2')
# f2 = FastqGeneralIterator(gzip.open(r2_fastq_file, 'rt')) 

barcodes_1 =  np.chararray(len(aln), itemsize=barcode1_length, unicode=True) 
# barcodes_2 =  np.chararray(len(aln), itemsize=barcode2_length, unicode=True) 
# The 0-position is a dummy. We use positions 1 and 2 for storing reads 1 and 2
reads = [None, None, None]
reads[1] = next(f1)
# reads[2] = next(f2)
# assert(reads[1][0] == reads[2][0])
read_count = 0
aln_indices = aln.index
for pos in range(len(aln_indices)):
    while ((int(reads[1][0]) != aln_indices[pos])):
        reads[1] = next(f1)
#         reads[2] = next(f2)
#         assert(reads[1][0] == reads[2][0])
        read_count += 1
        if read_count % 1e6 == 0:
            print(f'{read_count} reads processed.')
    barcodes_1[pos] = reads[barcode1_read][1][barcode1_start:barcode1_end]
#     barcodes_2[pos] = reads[barcode2_read][1][barcode2_start:barcode2_end]
print(f'Total {read_count} reads processed.')


# ## Count the number of reads for each barcode-insert pair and filter to barcode pairs that got at least 10 counts

aln['barcode_1'] = barcodes_1
# aln['barcode_2'] = barcodes_2


barcode_insert_counts = (aln
                         .groupby(['barcode_1', 'ref'])
                         .size()
                         .to_frame(name='read_count')
                         .reset_index()
                         .rename(columns={'ref': 'insert_num'})
                         .sort_values(by='read_count', ascending=False)
                         .reset_index(drop=True)
                         )
# barcode_insert_counts = barcode_insert_counts[barcode_insert_counts['read_count'] >= 10].reset_index(drop=True)
barcode_insert_counts.to_csv(output_file, sep = '\t', index_label = 'barcode_num')


(barcode_insert_counts
 .loc[barcode_insert_counts.read_count >= 4]
 .groupby('barcode_1')
 .size()
 .to_frame(name='n')
 .reset_index()
 .sort_values(by='n', ascending=False)
)


p = (
    pn.ggplot(barcode_insert_counts.reset_index(),
              pn.aes(x='index', y='read_count'))
    + pn.geom_line()
    + pn.labs(x="Barcode pair number", y="Read counts")
    + pn.scale_y_log10()
)

p.save(output_file.replace('.tsv.gz', '.png'), height=4, width=3)

p

