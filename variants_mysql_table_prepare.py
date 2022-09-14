#!/usr/bin/env python
# coding: utf-8

import os
import sys
import glob

try:
	infile = sys.argv[1]
except:
	sys.exit('<variants_table e.g Homo_sapiens.GRCh38.95.variants.intersections>')

outfile = '.'.join(infile.split('.')[:3]) + '.vrnts'

def rel_canvas_pos(elem_genomic_start, elem_genomic_end, utr_length):
	elem_genomic_start, elem_genomic_end, utr_length = int(elem_genomic_start), int(elem_genomic_end), int(utr_length)
	segment = round((1/float(utr_length) * 200),5)
	if elem_genomic_start == 1:
		#pos_start = round(((int(elem_genomic_start) * 200)/float(utr_length)),5)
		pos_start = 1
		if elem_genomic_end == utr_length:
			pos_end = 200.0
		else:
			pos_end = round((elem_genomic_end * segment),5)
		span = pos_end - pos_start
	elif elem_genomic_start == utr_length:
		pos_start = 200.0 
		pos_end = 0.0
		span = pos_end 
	elif elem_genomic_start == elem_genomic_end:
		pos_start = round(((int(elem_genomic_start) * 200)/float(utr_length)),5)
		pos_end = 0.0
		span = pos_end 
	else:
		pos_start = round(((int(elem_genomic_start) * 200)/float(utr_length)),5)
		if elem_genomic_end == utr_length:
			pos_end = 200 - pos_start
			span = pos_end
		else:
			pos_end = round((elem_genomic_end * segment),5)
			span = pos_end - pos_start
 
	#if str(pos_start).startswith('0.'): 
	#	pos_start = pos_start + 1 
	return (pos_start, span)

def genomic_pos(line):
	if line != []:
		chrom = line[0]
		utr_start = int(line[1])
		utr_end = int(line[2])
		variant_start = int(line[8])
		variant_end = int(line[9])
		overlap = int(line[-1])
		strand = line[5]
		utr_length = (utr_end - utr_start) + 1
		if line[5] == '+':
			if utr_start < variant_start and utr_end > variant_end: #phast interna
				genomic_variant_start = str((variant_start - utr_start)+1)
				genomic_variant_end = str(((variant_start - utr_start)+1) + overlap)
				#print 'variant_interna_' + str(line) + genomic_variant_start  + ' ' + genomic_variant_end + ' ' + str(variant_start), str(variant_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start < variant_start and utr_end < variant_end: #phast cavallo 3 primo
				genomic_variant_start = str((variant_start - utr_start)+1)
				genomic_variant_end = str(((variant_start - utr_start)+1) + overlap)
				#print 'variant_cavallo_3_primo_' + str(line), genomic_variant_start + ' ' + genomic_variant_end + ' ' + str(variant_start), str(utr_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(utr_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start > variant_start and utr_end > variant_end: #phast cavallo 5 primo
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap)
				#print 'variant_cavallo_5_primo_' + str(line) + genomic_variant_start + ' ' + genomic_variant_end + ' ' + str(utr_start), str(variant_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(utr_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start > variant_start and utr_end < variant_end: #phast esterna > utr_length
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap)
				#print 'variant_esterna' + str(line) + genomic_variant_start + ' ' + genomic_variant_end + ' ' + str(utr_start), str(utr_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(utr_start), str(utr_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start == variant_start and utr_end == variant_end:
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap)
				#print 'variant_uguale_utr_length' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(variant_start), str(variant_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start < variant_start and utr_end == variant_end:
				genomic_variant_start = str((variant_start - utr_start)+1)
				genomic_variant_end = str(((variant_start - utr_start)+1) + overlap)
				#print 'variant_uguale_al_3_primo_e_maggiore_del_5_primo' + str(line) +  genomic_variant_start  + ' ' + genomic_variant_end + ' ' + str(variant_start), str(variant_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start > variant_start and utr_end == variant_end:
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap) 
				#print 'variant_uguale_al_3_primo_e_minore_del_5_primo' + str(line) +  genomic_variant_start  + ' ' + genomic_variant_end + ' ' + str(utr_start), str(variant_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(utr_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start == variant_start and utr_end > variant_end:
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap) 
				#print 'variant_uguale_al_5_primo_e_minore_del_3_primo' + str(line) +  genomic_variant_start  + ' ' + genomic_variant_end + ' ' + str(variant_start), str(variant_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start == variant_start and utr_end < variant_end:
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap) 
				#print 'variant_uguale_al_5_primo_e_maggiore_del_3_primo' + str(line) +  genomic_variant_start  + ' ' + genomic_variant_end + ' ' + str(variant_start), str(utr_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(utr_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
		else:
			if utr_start < variant_start and utr_end > variant_end: #phast interna
				genomic_variant_start = str((utr_end - variant_end)+1)
				genomic_variant_end = str(((utr_end - variant_end)+1) + overlap)
				#print 'variant_interna_' + str(line) + genomic_variant_start + ' ' + genomic_variant_end + ' ' + str(variant_start), str(variant_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start > variant_start and utr_end > variant_end: #phast cavallo 3 primo
				genomic_variant_start = str((utr_end - variant_end)+1)
				genomic_variant_end = str(((utr_end - variant_end)+1) + overlap)
				#print 'variant_cavallo_3_primo_' + str(line), str((utr_end - variant_end)+1) + ' ' + str(((utr_end - variant_end)+1) + overlap) + ' ' + str(utr_start), str(variant_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(utr_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start < variant_start and utr_end < variant_end: #phast cavallo 5 primo
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap) 
				#print 'variant_cavallo_5_primo_' + str(line) + ' ' + str(1) + ' ' + str(1 + overlap) + ' ' + str(variant_start) + ' ' + str(utr_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(utr_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start > variant_start and utr_end < variant_end:
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap) 
				#print 'variant_esterna_' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(utr_start), str(utr_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(utr_start), str(utr_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start == variant_start and utr_end == variant_end:
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap) 
				#print 'variant_uguale_utr_length' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(variant_start), str(variant_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start == variant_start and utr_end < variant_end:
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap) 
				#print 'variant_uguale_al_3_primo_e_maggiore_del_5_primo' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(variant_start), str(utr_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(utr_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start == variant_start and utr_end > variant_end:
				genomic_variant_start = str((utr_end - variant_end)+1)
				genomic_variant_end = str(((utr_end - variant_end)+1) + overlap)
				#print 'variant_uguale_al_3_primo_e_minore_del_5_primo' + str(line) + str((utr_end - variant_end)+1) + ' ' + str(((utr_end - variant_end)+1) + overlap) + ' ' + str(utr_start), str(utr_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(utr_start), str(utr_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
			elif utr_start > variant_start and utr_end == variant_end:
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap) 
				return(genomic_variant_start, genomic_variant_end, chrom, str(utr_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))
				#print 'variant_uguale_al_5_primo_e_minore_del_3_primo' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(utr_start), str(variant_end)
			elif utr_start < variant_start and utr_end == variant_end:
				genomic_variant_start = str(1)
				genomic_variant_end = str(1 + overlap) 
				#print 'variant_uguale_al_5_primo_e_maggiore_del_3_primo' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(variant_start), str(variant_end)
				return(genomic_variant_start, genomic_variant_end, chrom, str(variant_start), str(variant_end), strand, rel_canvas_pos(genomic_variant_start, genomic_variant_end, utr_length))



dic_line = {} #1       ENSG00000142609 ENST00000642590 - [5_utr] [5_utr_variants] [3_utr] [3_utr_variants]       

with open(infile,'r') as f:
	for line in f:
		#line = line.strip('\n').split('\t')
		line = filter(None, (line.strip('\n').replace('\'','').split('\t')))   
		genomic_variants = genomic_pos(line)   
		dic_line.setdefault((line[0],line[3], line[4], line[5]), ([], [], [], []))
		if line[6] == '5_utr':
			dic_line.setdefault((line[0],line[3], line[4], line[5]))[0].append(line)
			dic_line.setdefault((line[0],line[3], line[4], line[5]))[1].append(genomic_variants)
		else:
			dic_line.setdefault((line[0],line[3], line[4], line[5]))[2].append(line)
			dic_line.setdefault((line[0],line[3], line[4], line[5]))[3].append(genomic_variants)


with open(outfile, 'w+') as g:
	g.write('{}\t{}\t{}\t{}\t{}\n'.format('Chrom','Gene_ID','Transcript_id','Strand','Variants'))
	for k,v in sorted(dic_line.items(), key= lambda x: int(x[0][0]) if x[0][0].isdigit() else x[0][0]):                  
		g.write('{}\t{}\n'.format('\t'.join(k), v)) 

cmd = 'rm -rf {}'
file_to_remove = glob.glob('*.tmp*')
print 'Removing temporary files..'
for elem in file_to_remove:
	os.system(cmd.format(elem))
print 'DONE!'
print '----------------------------' 
print 'Gzipping {}'.format(outfile)
cmd2 = 'gzip {}'.format(outfile)
os.system(cmd2)
print 'DONE!'
