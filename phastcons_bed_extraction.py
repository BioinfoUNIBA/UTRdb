#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import subprocess
import sys
import os

try:
	utrdb_tsv = sys.argv[1]
	intersecting_table = sys.argv[2]
except:
	sys.exit('<utrdb_tsv e.g HUMAN_hg38_5_UTR_seq_CAGES_seq.tsv, intersecting_phastcons_table e.g phastConsElements100way>')

df = pd.read_csv(utrdb_tsv, sep='\t')

temp_infl = '%s.bed' %(utrdb_tsv.strip('.tsv'))
temp_intersect = '%s.bed' %(intersecting_table.strip('.txt'))
phastconstype = intersecting_table.strip('.txt')


#output_file = '%s.phastConsElements' %(utrdb_tsv.strip('.tsv'))
output_file = '{}.{}'.format(utrdb_tsv.strip('.tsv'), phastconstype)


cmd = '''awk 'BEGIN{FS=OFS="\t"}{gsub("chr","",$2); print $2,$3+1,$4,$5,$6}' \
               %s > %s''' %(intersecting_table, temp_intersect)

os.system(cmd)

print '%s COMPLETE' %(intersecting_table)

with open(temp_infl,"w+") as out2:
	for i in range(0,len(df)): 
		if not pd.isnull(df.iloc[i,7]):
		    for utr_5 in eval(df.iloc[i,7]):
		        out2.writelines('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(utr_5[-1], utr_5[0], utr_5[1], \
		                                                                 df.iloc[i,1], df.iloc[i,4], \
		                                                                 utr_5[2], '5_utr'))#temp_infl.writelines
		if not pd.isnull(df.iloc[i,11]):
		    for utr_3 in eval(df.iloc[i,11]):
		         out2.writelines('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(utr_3[-1], utr_3[0], utr_3[1], \
		                                                                 df.iloc[i,1], df.iloc[i,4], \
		                                                                 utr_3[2], '3_utr'))#temp_infl.writelines

print '%s COMPLETE' %(utrdb_tsv)
print 'Intersecting tables...'
cmd2 = '''bedtools intersect -a {} -b {} -wo > {}.tmp'''.format(temp_infl, temp_intersect, output_file)
os.system(cmd2)
print 'Intersection COMPLETE!'
print 'Sorting...'
cmd3 = '''sort -k1,1 -k2,2n -V {}.tmp > {}.tmp2'''.format(output_file, output_file)
os.system(cmd3)
print 'Sorting COMPLETE! check {}.tmp2'.format(output_file) 
