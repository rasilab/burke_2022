#!/usr/bin/env python
# coding: utf-8

import os
import sys
import subprocess as sp
import re
import pandas as pd


spectrum_folder = '../data/spectrum_files/'

for input_file in os.listdir(spectrum_folder):
    input_path = os.path.join(spectrum_folder, input_file)
    output_path = os.path.join(spectrum_folder.replace('spectrum_files', 'sesca_deconv'), input_file)
    cmd = ['python', 'SESCA_deconv.py',
           '@spect', input_path,
           '@lib', '../libs/Map_BB_DS-dTSC3.dat',
           '@write', output_path,
           '@err', '2',
           '@rep', '100']
    sp.call(cmd)
    print(input_file, " processed.")


output_folder = '../data/sesca_deconv/'

table = dict()

for file in os.listdir(output_folder):
    if not file.endswith('.txt'):
        continue
    data = open(os.path.join(output_folder, file)).read().split('\n')
    sample = file.split('.')[0]
    alpha = float([re.search('[\.\d]+$', line).group(0) for line in data if '#        Alpha' in line][0])
    beta = float([re.search('[\.\d]+$', line).group(0) for line in data if '#         Beta' in line][0])
    coil = float([re.search('[\.\d]+$', line).group(0) for line in data if '#         Coil' in line][0])
    table[sample] = {'alpha': alpha, 'beta': beta, 'coil': coil}
    
table = (
    pd.DataFrame.from_dict(table, orient='index')
    .reset_index()
    .rename({'index': 'sample'}, axis=1)
    .sort_values(by = 'sample')
)

table.to_csv('../data/sesca_inferred_ss.csv', index=False)

table

