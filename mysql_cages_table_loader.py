#!/usr/bin/env python
import pandas as pd
from sqlalchemy import create_engine
#from sqlalchemy.dialects.mysql import LONGTEXT
import sys

try:
	input_file = sys.argv[1]
except:
	sys.exit('<input_tsv_gz_file>')

df = pd.read_csv(input_file, compression="gzip", sep="\t", low_memory=False)
#df.rename(columns = {'Gene_name/Gene_ID':'Gene_name_Gene_ID'}, inplace = True)
#df = df.sort_values(by=['Gene_name_Gene_ID'])

engine = create_engine("mysql+pymysql://{user}:{pw}@127.0.0.1:{port}/?charset=utf8".format(user="root",
                               pw="", port=3307))


with engine.connect() as connection:
	avail_db = connection.execute('SHOW DATABASES').fetchall()
	avail_db =  map(lambda x: str(x[0]), avail_db)
	if not 'Cages_tables' in avail_db:
		print 'Cages_tables db not available, I\'m gonna creating it..'
		connection.execute('CREATE DATABASE Cages_tables')

count_row = df.shape[0]
print 'Your df {} contains {} rows'.format(input_file, count_row)
print 'Going to upload its content...'

engine = create_engine("mysql+pymysql://{user}:{pw}@127.0.0.1:{port}/{db}?charset=utf8".format(user="root",
                               pw="", port=3307, db="Cages_tables"))

#dtype = {
#    "Editing_events": LONGTEXT
#}

sql_table = input_file.replace('.cages.gz','').replace('.','_').replace('-','_') + '_cages'
df.to_sql(sql_table, con=engine, if_exists='replace', index=False, chunksize=1000)#, dtype=dtype)
sql_count = 'SELECT COUNT(*) FROM {}'.format(sql_table).rstrip(',')
print 'Your uploaded utrs_table contains {} rows'.format(list(engine.execute(sql_count))[0][0])
print '{} Succesfully UPLOADED!'.format(sql_table) 
print 'Compressing {}'.format(sql_table)
cmd = "ALTER TABLE `{}` ROW_FORMAT=COMPRESSED;".format(sql_table)
engine.execute(cmd)
print 'Compression EXECUTED, ALL DONE for {}!'.format(sql_table)
print 'Indexing {}'.format(sql_table)
#cmd_index = "CREATE  UNIQUE INDEX transcript_id ON {} (Transcript_id(255))".format(sql_table)
cmd_index = "CREATE  INDEX transcript_id ON {} (Transcript_id(255))".format(sql_table)
engine.execute(cmd_index)
print 'Indexing EXECUTED, ALL DONE for {}!'.format(sql_table) 
engine.dispose()
print '-------------------------------------------------------------------------------------'
