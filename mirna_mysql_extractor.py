#!/usr/bin/env python
# coding: utf-8

# In[12]:


import MySQLdb
import pandas as pd
import itertools
import sys

try:
	input_tsv_file = sys.argv[1]
	release = sys.argv[2].strip()
	ensembl_db = sys.argv[3].strip()
except:
	sys.exit('<input_tsv_file, release, alldb|ensembl>') 

dic_url_port = {'alldb': ['mysql-eg-publicsql.ebi.ac.uk', 4157], 
				'ensembl': ['ensembldb.ensembl.org', 5306]}

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

def mirna_extractor(seq_coord): #, utr_type):
    mirna_targt = []
    dic_strand = {'-': '-1', '+': '1'}
    for cord in (eval(seq_coord)):
        seq_reg_start, seq_reg_end, seq_reg_strand, seq_reg_id = cord[0], cord[1],                                                                  cord[2], dic_chrom_id.get(cord[3], None) #dic_chrom_id[cord[3]]

        query = 'select accession, display_label, evidence, seq_region_id, seq_region_start, seq_region_end, seq_region_strand from mirna_target_feature where seq_region_id =\'%s\'         AND seq_region_start >=\'%s\' AND seq_region_end <=\'%s\' AND seq_region_strand =\'%s\';' %(seq_reg_id,
                                                                      seq_reg_start-20, 
                                                                      seq_reg_end+20, dic_strand[seq_reg_strand])

        cursor.execute(query)
        if any(cursor):
            mirna_targt.append(list(cursor.fetchall()))
        else:
            mirna_targt.append([])
    return mirna_targt

def overall_mirnas(seq_coord):
    mirna_targt = []
    dic_strand = {'-': '-1', '+': '1'}
    seq_coord = eval(seq_coord)
    seq_reg_start, seq_reg_end, seq_reg_strand, seq_reg_id = str(seq_coord[0][0]-20), str(seq_coord[-1][1]+20),                                                             seq_coord[-1][2], str(dic_chrom_id.get(seq_coord[-1][3], None))

    query = 'select accession, display_label, evidence, seq_region_id, seq_region_start, seq_region_end, seq_region_strand from mirna_target_feature where seq_region_id =\'%s\'     AND seq_region_start >=\'%s\' AND seq_region_end <=\'%s\' AND seq_region_strand =\'%s\';' %(seq_reg_id,
                                                              seq_reg_start, 
                                                              seq_reg_end, dic_strand[seq_reg_strand])
    cursor.execute(query)
    if any(cursor):
        mirna_targt = (list(cursor.fetchall()))
    return map(lambda x: tuple(str(i) for i in x), mirna_targt)

    
#db = MySQLdb.connect(host = 'ensembldb.ensembl.org', user = 'anonymous', port = 5306)
db = MySQLdb.connect(host= dic_url_port[ensembl_db][0], user = 'anonymous', port = dic_url_port[ensembl_db][1])
cursor = db.cursor()

core_db = cursor.execute('show databases LIKE "%s_core_%s%%";' %(input_tsv_file.lower().split('.')[0], release))
core_db = list(cursor.fetchall())[0][0]
funcgen_db = cursor.execute('show databases LIKE "%s_funcgen_%s%%";' %(input_tsv_file.lower().split('.')[0], release))
funcgen_db = list(cursor.fetchall())[0][0]

#sql = 'use homo_sapiens_core_95_38;'
#sql = 'use %s_core_95_38;' %(input_tsv_file.lower().split('.')[0])
sql = 'use %s;' %(core_db)
cursor.execute(sql)
df = pd.read_csv(input_tsv_file, sep = '\t')
chrom = list(set(df['Chrom']))
#chrom = [i for i in range(1,23)] + ['X','Y'] 
dic_chrom_id = {}
for chrm in chrom:
    cursor.execute('select seq_region_id from seq_region where name="%s" limit 1;' %(chrm))
    dic_chrom_id.setdefault(str(chrm), cursor.fetchone()[0])

#sql2 = 'use homo_sapiens_funcgen_95_38;'
#sql2 = 'use %s_funcgen_95_38;' %(input_tsv_file.lower().split('.')[0])
sql2 = 'use %s;' %(funcgen_db)
cursor.execute(sql2)
#df = pd.read_csv('HUMAN_hg38_5_UTR_seq_CAGES_seq.tsv', sep = '\t')
mirna_targets_table = results = cursor.execute("SELECT * from mirna_target_feature limit 1")
if mirna_targets_table:
	col_names = pd.Series(list(df.head())) #series from dataframe header
	drop_col_indexes = list(col_names[[5,6,8,9,10,12]])
	utrs_df = df.drop(drop_col_indexes, axis = 1)
	utrs_df.insert(6, '5_UTRs_exons_mirnas', '', allow_duplicates=False)
	utrs_df.insert(7, '5_UTRs_introns_mirnas', '', allow_duplicates=False)
	utrs_df.insert(9, '3_UTRs_exons_mirnas', '', allow_duplicates=False)
	utrs_df.insert(10, '3_UTRs_introns_mirnas', '', allow_duplicates=False)

	#print list(enumerate(utrs_df.columns))


	for i in range(0,len(utrs_df)):
		print 'Processing row %s genename: %s  transcript_id %s' %(i, utrs_df.iloc[i, 0], utrs_df.iloc[i, 4]) 
		utr5_coord = utrs_df.iloc[i, 5]
		utr3_coord = utrs_df.iloc[i, 8] #utrs_df.iloc[i, 7] 
		if type(utr5_coord) != float:
		    utr_5_mirnas = map(lambda x: [tuple(map(str, list(i))) for i in x], mirna_extractor(utr5_coord))
		    utrs_df.at[i, '5_UTRs_exons_mirnas'] = utr_5_mirnas
		    utrs_df.at[i, '5_UTRs_introns_mirnas'] = list(set(overall_mirnas(utr5_coord)).difference(set(list(itertools.chain(*utr_5_mirnas)))))
		if type(utr3_coord) != float:
		    utr_3_mirnas = map(lambda x: [tuple(map(str, list(i))) for i in x], mirna_extractor(utr3_coord))
		    utrs_df.at[i, '3_UTRs_exons_mirnas'] = utr_3_mirnas
		    utrs_df.at[i, '3_UTRs_introns_mirnas'] = list(set(overall_mirnas(utr3_coord)).difference(set(list(itertools.chain(*utr_3_mirnas)))))

	outfile = './' + input_tsv_file.replace('.utr_table.tsv','.mirnas')	    
	export_csv = utrs_df.to_csv(outfile, index = None, header=True, sep='\t') #Don't forget to add '.csv' at the end of the path
	print 'mirnas targets extraction COMPLETE!'
else:
	print 'mirnas targets not available for this organism!'
