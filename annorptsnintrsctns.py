#!/usr/bin/env python
import sys

try:
	annorepeats = sys.argv[1]
except:
	sys.exit('<.annorepeats>')

intersections = annorepeats.replace('annorepeats','annorptsnintrsctns')
with open(intersections, 'w+') as g, open(annorepeats,'r') as h:
        header = ['transcript_id', 'gene_id', 'chrom', 'strand', 'utr_start', 'utr_end', 'utr_type', 'repeat_intersect_start',
                'repeat_intersect_end', 'repeat_overlap_len', 'repeat_name', 'repeat_type','repeat_class']
        g.write('%s\n' %('\t'.join(header)))
        g.write('\n'.join(['\t'.join([c[4],c[3],c[0],c[5],c[1],c[2],c[6],c[8],c[9],c[14].strip(),c[13],c[10],c[11]])
                for c in [i.split('\t') for i in h.readlines()]]))
print 'JOB COMPLETED for %s!' %(annorepeats).upper()
