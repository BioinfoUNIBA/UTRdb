import sys, json, requests

try:
	infile = sys.argv[1]
	all_go = sys.argv[2]
except:
	sys.exit('<Homo_sapiens.GRCh38.95.gofnctns >')

outfile = infile.replace('gofnctns','gofnctns2')                            
def get_go(go):                                              
	requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query={}&limit=1&page=1".format(go)
	r = requests.get(requestURL, headers={ "Accept" : "application/json"})
	responseBody = r.text
	y = json.loads(responseBody)
	return(y['results'][0]['aspect'])


dic_go = {}                                           
with open(all_go,'r') as f:            
	for line in f:                                                                                             
		line = line.split('\t')                                           
		go = line[0]     
		descr = line[1].strip() 
		dic_go.setdefault(go, descr)

with open(infile,'r') as f, open(outfile, 'w') as g:
	new_head = '\t'.join([f.readline().strip(), 'GO_Aspect'])
	g.write('{}\n'.format(new_head))
	for line in f:
		line = line.strip().split('\t')
		go = line[2]
		if go != 'None':
			if not dic_go.get(go):
				line = line + [(get_go(go))]
			else:
				line = line + [dic_go.get(go)]
			g.write('{}\n'.format('\t'.join(line)))
		else:
			line = line + ['None']
			g.write('{}\n'.format('\t'.join(line)))

