#!/usr/bin/env python2.7
# coding: utf-8

import pandas as pd
import subprocess
import sys
import os

try:
	utrdb_tsv = sys.argv[1]
	intersecting_table = sys.argv[2]
except:
	sys.exit('<utrdb_tsv e.g HUMAN_hg38_5_UTR_seq_CAGES_seq.tsv, intersecting_variants_table e.g Homo_sapiens.GRCh38.95.variants>')

df = pd.read_csv(utrdb_tsv, sep='\t')

temp_infl = '%s.bed' %(utrdb_tsv.strip('.tsv'))
temp_intersect = '%s.bed' %(intersecting_table)
temp_intersection_file = temp_intersect.replace('.bed','.intersections.bed')

with open(intersecting_table,'r') as d:
	if not d.readlines()[2].split('\t')[0].isdigit(): #if chr1 or 1 apply corresponding awk command
		cmd = '''awk 'BEGIN{FS=OFS="\t"}{gsub("chr","",$1); if (NR!=1) print $1,$2,$2,$3,$4,$5,$6,$7,$8,$9,$12}' \
						%s > %s''' %(intersecting_table, temp_intersect)
	else:
		cmd = '''awk 'BEGIN{FS=OFS="\t"}{if (NR!=1) print $1,$2,$2,$3,$4,$5,$6,$7,$8,$9,$12}' \
						%s > %s''' %(intersecting_table, temp_intersect)

print 'PREPARING %s from %s' %(temp_intersect, intersecting_table)
os.system(cmd)
print '%s COMPLETE' %(temp_intersect)

print 'Removing duplicate rows from {}'.format(temp_intersect)
cmd2 =  '''awk '!visited[$0]++' {} > {}'''.format(temp_intersect, temp_intersect.replace('.bed','.dedupl.bed'))
os.system(cmd2)
print '%s COMPLETE' %(temp_intersect.replace('.bed','.dedupl.bed'))

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

print '%s COMPLETE' %(temp_infl)

print 'Intersecting tables...'
cmd3 = '''bedtools intersect -a {} -b {} -wo > {}.tmp'''.format(temp_infl, temp_intersect.replace('.bed','.dedupl.bed'), temp_intersection_file.replace('.bed','.dedupl.bed'))
os.system(cmd3)
print 'Intersection COMPLETE!'

print 'Sorting...'
cmd4 = '''sort -k1,1 -k2,2n -V {}.tmp > {}'''.format(temp_intersection_file.replace('.bed','.dedupl.bed'), temp_intersection_file.replace('.bed',''))
os.system(cmd4)
print 'Sorting COMPLETE! check {}'.format(temp_intersection_file.replace('.bed',''))
