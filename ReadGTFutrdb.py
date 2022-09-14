#!/usr/bin/env python2.7
import sys,os

try:
	infile=sys.argv[1] #'gene.gtf' #'gencode.v33lift37.annotation.gtf'
except: 
	sys.exit('<gtf file>')


org=infile.split('.')[1]
tmpfile = infile.replace('.gtf','.tmp')
outfile = tmpfile.replace('.tmp','.annotations')

#edifile='listEditingSites.txt'

def getL(last):
	d={}
	last=last.replace('"','')
	last=last.rstrip(';')
	for i in last.split(';'):
		if i.lstrip()=='': continue
		ll=(i.lstrip()).split()
		d[ll[0].replace('bio','')]=ll[1]
	return d
	
def makeCluster(allcoord):
	cluster=[]
	remaining=[]
	c1=allcoord[0][0]
	c2=allcoord[0][1]
	for i in range(len(allcoord)):
		if allcoord[i]!=(c1,c2):
			if c1<=allcoord[i][0]<=c2:
				cluster.append(allcoord[i])
				if allcoord[i][1]>c2:
					c2=allcoord[i][1]
			else:
				remaining.append(allcoord[i])
		else:
			cluster.append((c1,c2))
	return (c1,c2),remaining

def Clusters(interval):
	coords=[]
	interval.sort()
	while len(interval)!=0:
		coord,interval=makeCluster(interval)
		coords.append(coord)
	return coords
	
d={}
dg={}
dex={}
f=open(infile)
for i in f:
	if i.startswith('#'): continue
	l=(i.strip()).split('\t')
	print l
	if l[2] not in ['exon','CDS']: continue
	last=getL(l[-1])
	if (last['gene_id'].count('.') == 1 and last['gene_id'].split('.')[-1].isdigit()):
		gene=last['gene_id'].split('.')[0]
	else:
		gene=last['gene_id']
	if (last['transcript_id'].count('.') == 1 and last['transcript_id'].split('.')[-1].isdigit()):
		tr=last['transcript_id'].split('.')[0]
	else:
		tr=last['transcript_id']
	try: gname=last['gene_name']
	except: gname=gene
	try: gtype=last['gene_type']
	except: gtype='unknown'
	try: ttype=last['transcript_type']
	except: ttype='unknown'
	dg[gene]=(gname,gtype)
	if d.has_key(gene):
		if d[gene].has_key(tr):
			d[gene][tr].append((int(l[3]),int(l[4]),l[2]))
		else:
			d[gene][tr]=[(int(l[3]),int(l[4]),l[2])]
			d[gene]['info'][2][tr]=ttype
	else:
		d[gene]={}
		d[gene]['info']=[l[0],l[6],{tr:ttype}]
		d[gene][tr]=[(int(l[3]),int(l[4]),l[2])]
	if dex.has_key(gene):
		if l[2]=='exon': dex[gene].append((int(l[3]),int(l[4])))
	else:
		if l[2]=='exon': dex[gene]=[(int(l[3]),int(l[4]))]
f.close()

lg=[]
with open(tmpfile, 'w') as g:
	header =  ['geneid', 'symbol', 'gtype', 'trid', 'ttype', 'startg', 'endg', 'start', 'end', 'strand', 'chr', 'feat', 'org']  
	g.write('{}\n'.format('\t'.join(header)))                     
	for i in d:
		#print i,dg[i]
		allex=dex[i]
		allex.sort()
		gcoord=(allex[0][0],allex[-1][1])
		strand=d[i]['info'][1]
		chr=d[i]['info'][0]
		for j in d[i].keys():
			if j=='info': continue
			trt=d[i]['info'][2][j]
			cds=[x for x in d[i][j] if x[2]=='CDS']
			cds.sort()
			if len(cds)>0: cdsc=(cds[0][0],cds[-1][1])
			else: cdsc=()
			str_tr=[]
			if len(cdsc)>0:
				# gene ha cds
				up,md,dw=[],[],[]
				stex=[]
				for k in d[i][j]:
					if k[2]!='exon': continue
					if cdsc[0] > k[1] and cdsc[1] > k[1]: up.append(k)
					if cdsc[1] < k[0] and cdsc[1] < k[1]: dw.append(k)
					if cdsc[0]< k[0]< cdsc[1] and cdsc[0]< k[1]< cdsc[1]: md.append(k)
					if k[0]<=cdsc[0]<=k[1]: stex.append(k)
					if k[0]<=cdsc[1]<=k[1]: stex.append(k)
				stex.sort()
				if stex[0][0] == cdsc[0]: 
					md.append(stex[0])
				elif stex[0][0] < cdsc[0] < stex[0][1]:
					up.append((stex[0][0],cdsc[0]-1,''))
					md.append((cdsc[0],stex[0][1],''))
				elif stex[0][1] == cdsc[0]: up.append(stex[0])
				if len(stex) > 1:
					if stex[1][0] == cdsc[1]: 
						dw.append(stex[1])
					elif stex[1][0] < cdsc[1] < stex[1][1]:
						md.append((stex[1][0],cdsc[1],''))
						dw.append((cdsc[1]+1,stex[1][1],''))
					elif stex[1][1] == cdsc[1]: 
						md.append(stex[1])
				if strand=='+':
					for k in up: str_tr.append((k[0],k[1],'5UTR',strand,'red'))
					for k in md: str_tr.append((k[0],k[1],'CDS',strand,'orange'))
					for k in dw: str_tr.append((k[0],k[1],'3UTR',strand,'blue'))
				elif strand=='-':
					for k in dw: str_tr.append((k[0],k[1],'5UTR',strand,'red'))
					for k in md: str_tr.append((k[0],k[1],'CDS',strand,'orange'))
					for k in up: str_tr.append((k[0],k[1],'3UTR',strand,'blue'))
				str_tr.sort()
				for k in range(len(str_tr)-1):
					c1=str_tr[k][1]+1
					c2=str_tr[k+1][0]-1
					if (c2-c1)+1 > 0: 
						str_tr.append((c1,c2,'Intron',strand,'black'))
				str_tr.sort()
				#print str_tr 
				for k in str_tr:
					line=[i,dg[i][0],dg[i][1],j,trt,str(gcoord[0]),str(gcoord[1]),str(k[0]),str(k[1]),k[3],chr,k[2],org]
					#line="{x: %i,x2: %i, y: 0, pointWidth: 20,color: '%s',name:'Gene: %s - %s<br>Transcript: %s - %s<br>Strand: %s<br>Feature: %s<br>%s'}" %(k[0],k[1],k[4],dg[i][0],dg[i][1],j,trt,k[3],k[2],chr)
					#print dg[i][0]+'\t'+j+'\t'+line
					#print '\t'.join(line)
					g.write('{}\n'.format('\t'.join(line)))
			else:
				color='green'
				feat='EXON'
				exon=[x for x in d[i][j] if x[2]=='exon']
				exon.sort()
				intron=[(exon[x][1]+1,exon[x+1][0]-1) for x in range(len(exon)-1)]
				for k in exon:
					line=[i,dg[i][0],dg[i][1],j,trt,str(gcoord[0]),str(gcoord[1]),str(k[0]),str(k[1]),strand,chr,feat,org]
					#print '\t'.join(line)
					g.write('{}\n'.format('\t'.join(line)))
				for k in intron:
					line=[i,dg[i][0],dg[i][1],j,trt,str(gcoord[0]),str(gcoord[1]),str(k[0]),str(k[1]),strand,chr,'Intron',org]
					#line="{x: %i,x2: %i, y: 0, pointWidth: 20,color: '%s',name:'Gene: %s - %s<br>Transcript: %s - %s<br>Strand: %s<br>Feature: %s<br>%s'}" %(k[0],k[1],color,dg[i][0],dg[i][1],j,trt,strand,feat,chr)
					#print dg[i][0]+'\t'+j+'\t'+line
					#print '\t'.join(line)
					g.write('{}\n'.format('\t'.join(line)))


cmd = '''awk '!visited[$0]++' {} > {}'''.format(tmpfile, outfile)
print 'Removing duplicates...'
os.system(cmd)
print 'Removing temporary files'
cmd1 = '''rm -rf {}'''.format(tmpfile)
os.system(cmd1)
print 'Gzipping annotations file'
cmd2 = '''gzip {}'''.format(outfile)
os.system(cmd2)

