#!/usr/bin/env python

import os
import pandas as pd
import subprocess
import tempfile
import sys

try:
    utrdb_tsv = sys.argv[1]
    intersecting_table = sys.argv[2]
except:
    sys.exit('<utrdb_tsv e.g HUMAN_hg38_5_UTR_seq_CAGES_seq.tsv, intersecting_table e.g rediportal_hg19liftedhg38>')

df = pd.read_csv(utrdb_tsv, sep='\t')

temp_infl = '%s.bed' %(utrdb_tsv.strip('.tsv'))
temp_intersect = '%s' %(intersecting_table.replace('.txt', '.bed'))
output_file = output_file = '{}.{}'.format(utrdb_tsv.strip('.tsv'), 'rdprtl')

#print temp_intersect

cmd = '''awk 'BEGIN{FS=OFS="\t"}{gsub("chr","",$1); print $1,$2,$3,$5,$6,$4,$7}' %s > %s ''' %(intersecting_table, temp_intersect)
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
cmd2 = '''awk '!visited[$0]++' %s > %s''' %(temp_infl, temp_infl.replace('bed','dedup.bed'))
os.system(cmd2)
#print cmd2
print 'Removing duplicates from {}'.format(temp_infl)

print 'Intersecting tables...'
cmd3 = '''bedtools intersect -a {} -b {} -wo > {}.tmp'''.format(temp_infl.replace('bed','dedup.bed'), temp_intersect, output_file)
os.system(cmd3)
#print cmd3

print 'INTERSECTION COMPLETE!'

print 'Sorting...'
cmd4 = '''sort -k1,1 -k2,2n -V {}.tmp > {}.tmp2'''.format(output_file, output_file)
os.system(cmd4)
#print cmd4
print 'Sorting COMPLETE! check {}.tmp2'.format(output_file) 
