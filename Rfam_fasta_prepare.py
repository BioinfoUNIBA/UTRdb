#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys
import os

try:
	infile = sys.argv[1]
except:
	sys.exit('<infile eg. HUMAN_hg38_5_UTR_seq_CAGES_seq.tsv>')

outfile = os.getcwd() + '/' + infile.replace('utr_table.tsv','fa').split('/')[-1]
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

df = pd.read_csv(infile, sep = '\t')
#col_names = pd.Series(list(df.head())) #series from dataframe header

with open(outfile,'w') as f:
	for i in range(0,len(df)):
		print 'Processing row %s gene %s transcript %s' %(i, df.iloc[i,0], df.iloc[i,4])
		if type(df.iloc[i,7]) != float:
		    for utr5_cord, seq5utr in zip(eval(str(df.iloc[i,7])),eval(str(df.iloc[i,8]))):
		        f.write('%s\n' %('>' + '5' + '_' + (df.iloc[i,4]) + '_' + str(utr5_cord).replace(', ','_').replace(')','').replace('(','')))
		        f.write('%s\n' %(seq5utr))
		if type(df.iloc[i,11]) != float:
		    for utr3_cord, seq3utr in zip(eval(str(df.iloc[i,11])),eval(str(df.iloc[i,12]))):
		        f.write('%s\n'  %('>' + '3' + '_' + (df.iloc[i,4]) + '_' + str(utr3_cord).replace(', ','_').replace(')','').replace('(','')))
		        f.write('%s\n' %(seq3utr))

print '%s CREATED!' %(outfile)

