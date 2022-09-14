#!/usr/bin/env python

import sys
import os
import cgi

try:
	cmscan = sys.argv[1]
	tblout = sys.argv[2]
except:
	sys.exit('<cmscan, tblout eg. \'Homo_sapiens.GRCh38.95.cmscan\', \'Homo_sapiens.GRCh38.95.rfam.tblout\'>')

outfile = '.'.join(tblout.split('.')[:-1])

def transcript_check(lista):
	a = lista[1:]
	b = map(lambda x: x.isdigit(), a).index(True)+1
	return(lista[1:b], lista[b:])



with open(cmscan, 'r') as e:
	content = ''.join(e.readlines()).split('//\n') 


dic_utr_cmscan = {}                                       
for i in content:                                                                                                                
	if i.find('>>') != -1:                                                                                                 
		i = '\n'.join(i.split('\n')[:-2])                            
		dic_utr_cmscan.setdefault(i.split('\n')[0].split()[1].replace('\'',''), []).append('<html><body><pre>' + cgi.escape(i) + '</pre></body></html>')


with open(tblout, 'r') as f, open(outfile, 'w') as g:
	header = ['Chrom', 'Transcript_id', 'utr_type', 'genomic_start', 'genomic_end', 'genomic_strand', 'rfam_genomic_start', 'rfam_genomic_end', 
				'idx', 'target_name', 'target_accession', 'query_name', 'query_accession', 
				'clan_name', 'mdl', 'mdl_from',  'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc', 'pass', 'gc', 'bias', 
				'score', 'E_value', 'inc', 'olp', 'anyidx', 'afrct1', 'afrct2', 'winidx', 'wfrct1', 'wfrct2', 'description_of_target', 'cmscan_output']
	g.write('{}\n'.format('\t'.join(header)))
	for line in f:
		if not line.startswith('#'):
			line = line.replace('\'','').strip().split()
			utr = line[3].split('_')
			utr_type = utr[0] + '_utr'
			transcript_id = '_'.join(transcript_check(utr)[0])
			
			genomic_start = transcript_check(utr)[1][0] #utr[2]
			genomic_end = transcript_check(utr)[1][1] #utr[3]
			genomic_strand = transcript_check(utr)[1][2] #utr[4].strip('"')
			chrom = transcript_check(utr)[1][3] #utr[5]
			if genomic_strand == '+':
				rfam_genomic_start = int(genomic_start) + (int(line[9]) - 1)
				rfam_genomic_end = int(genomic_start) + (int(line[10]) - 1)
			else:
				rfam_genomic_start = int(genomic_end) - int(line[10]) + 1
				rfam_genomic_end = int(genomic_end) - int(line[9]) + 1
			g.write('{}\n'.format('\t'.join([chrom, transcript_id, utr_type, genomic_start, genomic_end, genomic_strand, str(rfam_genomic_start), str(rfam_genomic_end)] + line[:26] + \
				[' '.join(line[26::])]) + '\t' + str(dic_utr_cmscan.get(line[3]))))


print 'Gzipping {}'.format(outfile)
cmd2 = 'gzip {}'.format(outfile)
os.system(cmd2)
print 'DONE!'	

