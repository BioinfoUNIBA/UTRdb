import sys
from lxml import etree   
parser = etree.XMLParser(recover=True)    

try:
	ensemble_ids = sys.argv[1]
	go_html = sys.argv[2]
except:
	sys.exit('<ensembl_ids, go_html')
     
def extractor(output):
	root = etree.fromstring(output, parser=parser)
	output = []
	for node in root:
		go_id, go_term = node.get('display_id'), node.get('description')
		line = '{}\t{}'.format(go_id, go_term)
		output.append(line)
	return list(set(output))

header = ['Gene_id','Transcript_id','GO_ID', 'GO_Term']
with open(ensemble_ids,'r') as a, open(go_html,'r') as b:
	ids = a.readlines()
	html = ''.join(b.readlines()).split('//')
print '\t'.join(header)
for i,j in zip(ids, filter(lambda x: x != '', map(str.strip, html))):
	if j.find('exceeded') == -1:
		if extractor(j[j.find('<opt>'):]):
			for go in extractor(j):
				print '{}\t{}'.format(i.strip(),go) 

#	for i,j in zip(ids,html):
#		print i,j
#	for i,j in zip(ids,html):
#		print i.strip(), ''.join(extractor(j))

