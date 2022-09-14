import pandas as pd
import sys
import os
import commands

try:
		infile = sys.argv[1]
except:
		sys.exit('<utr_table_tsv eg. Homo_sapiens.GRCh38.95.utr_table.tsv>')

df = pd.read_csv(infile, sep='\t')
a = df[['Gene_ID', 'Transcript_id']]

ensemble_ids = infile.replace('utr_table.tsv','ensemble_ids')
go_html = infile.replace('utr_table.tsv','go_html')

for i in range(0, len(a)):
	gene_id = a.iloc[i,0]
	transcript_id = a.iloc[i,1]
	print('Processing %s line %s/%s' %(transcript_id, i, len(a)))
	server = "https://rest.ensembl.org"
	ext = "/xrefs/id/%s?external_db=GO;all_levels=1;" %(transcript_id)
	cmd = "echo %s$'\t'%s >> %s" %(gene_id, transcript_id, ensemble_ids)
	os.system(cmd)
	i = 0
	while i == 0:
		if commands.getoutput('curl -o /dev/null -s -H "Accept: text/xml" -w \"%{http_code}\" \"' + server+ext +'"') in ['200','400']:
			i=1
	commands.getoutput('curl -Ss -H "Accept: text/xml" --connect-timeout 900 --retry 10 --retry-delay 15 "%s" >> %s' %(server+ext, go_html))
	cmd2 = 'echo // >> %s' %(go_html)
	os.system(cmd2)
print 'COMPLETE!'



