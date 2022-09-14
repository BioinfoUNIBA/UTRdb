import commands#requests
import pandas as pd
import sys
import time

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

try:
	infile = sys.argv[1]
except:
	sys.exit('<Utr_table.tsv>')

df = pd.read_csv(infile, sep = '\t')
df = df.drop_duplicates(subset ="Gene_ID", 
					 keep = "first", inplace = False) 
total_rows = df.shape[0]
logfile = infile.rstrip('.tsv') + '.logfile'
outfile = infile.rstrip('.tsv') + '.orthologs'
with open(logfile,'a', 0) as f, open(outfile,'a', 0) as g:
	g.write('%s\t%s\t%s\t%s\t%s\n' %('Gene_id', 'Symbol','Taxonomy', 'Species', 'Orthologue'))
	for row in range(0,len(df)):
			gene_name_id = df.iloc[row,0]
			gene_id = df.iloc[row,1]
			server = "https://rest.ensembl.org"
			ext = "/homology/id/%s?type=orthologues;format=condensed;content-type=text/xml" %(gene_id)
			f.write("Processing %s %.3f%%\n" %(gene_id, (float(row)/float(total_rows))*100))
			i = 0
			while i == 0:
				if commands.getoutput('curl -o /dev/null -s -H "Accept: text/xml" -w \"%{http_code}\" \"' + server+ext +'"') in ['200','400']:
					i=1
			output=commands.getoutput('curl -Ss -H "Accept: text/xml" --connect-timeout 900 --retry 10 --retry-delay 15 "%s"' %(server+ext))
			for i in (map(str, output.split('\n'))):
				i = i.strip().strip('<').strip('/>')
				if i.startswith('homologies'): 
					i = i.split() 
					if len(i) == 7:												
						species = i[4].split('=')[-1].strip('"')		  
						taxonomy_level = i[5].split('=')[-1].strip('"')  
						ortholog = i[1].split('=')[-1].strip('"')		 
						g.write('%s\t%s\t%s\t%s\t%s\n' %
						  			(gene_id, gene_name_id, taxonomy_level, species, ortholog))
			#time.sleep(5)  
	f.write('Orthologs extraction COMPLETE!\n')		




