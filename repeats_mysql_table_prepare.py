#!/usr/bin/env python
# coding: utf-8

import os
import sys
import glob

try:
	infile = sys.argv[1]
except:
	sys.exit('<sorted_annorepeats_table e.g Homo_sapiens.GRCh38.95.sorted.annorepeats>')

outfile = '.'.join(infile.split('.')[:-2]) + '.rpts'

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
		chrom = line[2]
		utr_start = int(line[4])
		utr_end = int(line[5])
		elem_start = int(line[7])
		elem_end = int(line[8])
		overlap = int(line[9])
		elem_name, elem_type, elem_class = line[10:]
		strand = line[3]
		utr_length = (utr_end - utr_start) + 1
		if line[3] == '+':
			if utr_start < elem_start and utr_end > elem_end: #phast interna
				genomic_elem_start = str((elem_start - utr_start)+1)
				genomic_elem_end = str(((elem_start - utr_start)+1) + overlap)
				#print 'elem_interna_' + str(line) + genomic_elem_start  + ' ' + genomic_elem_end + ' ' + str(elem_start), str(elem_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start < elem_start and utr_end < elem_end: #phast cavallo 3 primo
				genomic_elem_start = str((elem_start - utr_start)+1)
				genomic_elem_end = str(((elem_start - utr_start)+1) + overlap)
				#print 'elem_cavallo_3_primo_' + str(line), genomic_elem_start + ' ' + genomic_elem_end + ' ' + str(elem_start), str(utr_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(utr_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start > elem_start and utr_end > elem_end: #phast cavallo 5 primo
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap)
				#print 'elem_cavallo_5_primo_' + str(line) + genomic_elem_start + ' ' + genomic_elem_end + ' ' + str(utr_start), str(elem_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(utr_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start > elem_start and utr_end < elem_end: #phast esterna > utr_length
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap)
				#print 'elem_esterna' + str(line) + genomic_elem_start + ' ' + genomic_elem_end + ' ' + str(utr_start), str(utr_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(utr_start), str(utr_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start == elem_start and utr_end == elem_end:
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap)
				#print 'elem_uguale_utr_length' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(elem_start), str(elem_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start < elem_start and utr_end == elem_end:
				genomic_elem_start = str((elem_start - utr_start)+1)
				genomic_elem_end = str(((elem_start - utr_start)+1) + overlap)
				#print 'elem_uguale_al_3_primo_e_maggiore_del_5_primo' + str(line) +  genomic_elem_start  + ' ' + genomic_elem_end + ' ' + str(elem_start), str(elem_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start > elem_start and utr_end == elem_end:
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap) 
				#print 'elem_uguale_al_3_primo_e_minore_del_5_primo' + str(line) +  genomic_elem_start  + ' ' + genomic_elem_end + ' ' + str(utr_start), str(elem_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(utr_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start == elem_start and utr_end > elem_end:
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap) 
				#print 'elem_uguale_al_5_primo_e_minore_del_3_primo' + str(line) +  genomic_elem_start  + ' ' + genomic_elem_end + ' ' + str(elem_start), str(elem_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start == elem_start and utr_end < elem_end:
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap) 
				#print 'elem_uguale_al_5_primo_e_maggiore_del_3_primo' + str(line) +  genomic_elem_start  + ' ' + genomic_elem_end + ' ' + str(elem_start), str(utr_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(utr_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
		else:
			if utr_start < elem_start and utr_end > elem_end: #phast interna
				genomic_elem_start = str((utr_end - elem_end)+1)
				genomic_elem_end = str(((utr_end - elem_end)+1) + overlap)
				#print 'elem_interna_' + str(line) + genomic_elem_start + ' ' + genomic_elem_end + ' ' + str(elem_start), str(elem_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start > elem_start and utr_end > elem_end: #phast cavallo 3 primo
				genomic_elem_start = str((utr_end - elem_end)+1)
				genomic_elem_end = str(((utr_end - elem_end)+1) + overlap)
				#print 'elem_cavallo_3_primo_' + str(line), str((utr_end - elem_end)+1) + ' ' + str(((utr_end - elem_end)+1) + overlap) + ' ' + str(utr_start), str(elem_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(utr_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start < elem_start and utr_end < elem_end: #phast cavallo 5 primo
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap) 
				#print 'elem_cavallo_5_primo_' + str(line) + ' ' + str(1) + ' ' + str(1 + overlap) + ' ' + str(elem_start) + ' ' + str(utr_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(utr_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start > elem_start and utr_end < elem_end:
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap) 
				#print 'elem_esterna_' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(utr_start), str(utr_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(utr_start), str(utr_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start == elem_start and utr_end == elem_end:
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap) 
				#print 'elem_uguale_utr_length' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(elem_start), str(elem_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start == elem_start and utr_end < elem_end:
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap) 
				#print 'elem_uguale_al_3_primo_e_maggiore_del_5_primo' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(elem_start), str(utr_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(utr_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start == elem_start and utr_end > elem_end:
				genomic_elem_start = str((utr_end - elem_end)+1)
				genomic_elem_end = str(((utr_end - elem_end)+1) + overlap)
				#print 'elem_uguale_al_3_primo_e_minore_del_5_primo' + str(line) + str((utr_end - elem_end)+1) + ' ' + str(((utr_end - elem_end)+1) + overlap) + ' ' + str(utr_start), str(utr_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(utr_start), str(utr_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
			elif utr_start > elem_start and utr_end == elem_end:
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap) 
				return(genomic_elem_start, genomic_elem_end, chrom, str(utr_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))
				#print 'elem_uguale_al_5_primo_e_minore_del_3_primo' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(utr_start), str(elem_end)
			elif utr_start < elem_start and utr_end == elem_end:
				genomic_elem_start = str(1)
				genomic_elem_end = str(1 + overlap) 
				#print 'elem_uguale_al_5_primo_e_maggiore_del_3_primo' + str(line) + str(1) + ' ' + str(1 + overlap) + ' ' + str(elem_start), str(elem_end)
				return(genomic_elem_start, genomic_elem_end, chrom, str(elem_start), str(elem_end), strand, elem_name, elem_type, elem_class, rel_canvas_pos(genomic_elem_start, genomic_elem_end, utr_length))



dic_line = {} #1       ENSG00000142609 ENST00000642590 - [5_utr] [5_utr_phastcons] [3_utr] [3_utr_phastcons]       
with open(infile,'r') as f:
	for line in f:
		if line.startswith('transcript_id'):
			continue
		line = line.strip('\n').split('\t')   
		line_2 = [line[2], line[4], line[5], line[1], line[0], line[3], line[6], line[2], line[7], line[8], line[9], line[10], line[11], line[12]] #riordino la line
		genomic_repeats = genomic_pos(line)  
		dic_line.setdefault((line[2],line[1], line[0], line[3]), ([], [], [], [])) 
		if line[6] == '5_utr':
			dic_line.setdefault((line[2],line[1], line[0], line[3]))[0].append(line_2)
			dic_line.setdefault((line[2],line[1], line[0], line[3]))[1].append(genomic_repeats)
		else:
			dic_line.setdefault((line[2],line[1], line[0], line[3]))[2].append(line_2)
			dic_line.setdefault((line[2],line[1], line[0], line[3]))[3].append(genomic_repeats)


with open(outfile, 'w+') as g:
	g.write('{}\t{}\t{}\t{}\t{}\n'.format('Chrom','Gene_ID','Transcript_id','Strand','Rpts_name_type_class'))
	for k,v in sorted(dic_line.items(), key= lambda x: int(x[0][0]) if x[0][0].isdigit() else x[0][0]):                  
		g.write('{}\t{}\n'.format('\t'.join(k), v)) 

print 'Gzipping {}'.format(outfile)
cmd2 = 'gzip {}'.format(outfile)
os.system(cmd2)
print 'DONE!'
