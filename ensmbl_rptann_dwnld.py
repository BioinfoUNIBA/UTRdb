#! /usr/bin/env python
import os
import subprocess
import sys

try:
	infile = sys.argv[1]
	release = sys.argv[2]
	dtbs = int(sys.argv[3])
except:
	sys.exit('<infile_tsv, release, db[0->Ensembl|1->Ensembl Genomes all databases]>')

if dtbs == 0:
	database = "ensembldb.ensembl.org" 
	port = 5306
else:
	database = "mysql-eg-publicsql.ebi.ac.uk"
	port = 4157

infile_name = infile.split('.')[0].lower()

db = subprocess.check_output("mysql -N --host=%s --port=%s --user=anonymous -e \
	\'show databases WHERE `Database` LIKE \"%s_core_%s%%\";\'" \
	%(database, port, infile_name, release), shell=True).strip().split()[0]

print 'Using db %s' %(db)
outfile = infile.strip('.utr_table.tsv')

os.system("mysql -N --host=%s --port=%s --user=anonymous -e \
	'use '%s'; SELECT DISTINCT seq_region.name, repeat_feature.seq_region_start, \
	repeat_feature.seq_region_end, repeat_feature.seq_region_strand, repeat_consensus.repeat_type, \
	repeat_consensus.repeat_class, repeat_consensus.repeat_name FROM repeat_feature \
	JOIN repeat_consensus ON repeat_feature.repeat_consensus_id = repeat_consensus.repeat_consensus_id \
	JOIN seq_region ON repeat_feature.seq_region_id = seq_region.seq_region_id where \
	repeat_feature.seq_region_strand != 0;' > '%s.repeats'" \
	%(database, port, db, outfile))

print 'EXTRACTION COMPLETE!'
