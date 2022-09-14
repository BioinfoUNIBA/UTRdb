import pysam
import sys

try:
     infile = sys.argv[1]
     fasta_file = sys.argv[2]
except:
     sys.exit('<gtf_infile, fasta_file>')


pysam_object = pysam.FastaFile(fasta_file)

def dic_gtf_column_val(gtf_line):
    gtf_line = gtf_line.replace('"', '').rstrip(';')
    attr_value = {}
    for column in gtf_line.split(';'):
        attr = (column.lstrip()).split()  # split any attr internally by space
        attr_value.setdefault(attr[0], ' '.join(attr[1:]))  # voice eg gene_id, value
    return attr_value

def get_introns(exons_list):
    introns = []
    coords_prev_exons = map(lambda x: x[0], exons_list)
    coords_next_exons = map(lambda x: x[1], exons_list)
    chrom_strand_numb = map(lambda x: x[2:], exons_list)
    strand = map(lambda x: x[2], exons_list)
    if len(set(strand)) == 1:
        if strand[0] == '+':
            for x, y, z in zip(coords_next_exons, coords_prev_exons[1:], chrom_strand_numb):
                introns.append((x + 1, y - 1, z[0], z[1], z[2]))
        else:
            for x, y, z in zip(coords_next_exons[1:], coords_prev_exons, chrom_strand_numb):
                introns.append((x + 1, y - 1, z[0], z[1], z[2]))
    else:
        introns = []
    return introns

def get_sequences(coords):
    sequences = []
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for coord in coords:
        start = int(coord[0])
        end = int(coord[1])
        strand = coord[2]
        chrom = coord[3]
        seq = pysam_object.fetch(chrom, start - 1, end)
        if strand != '+':
            seq = "".join(complement.get(base, base) for base in reversed(seq))
        sequences.append(seq)
    return sequences

def swap_coords(exons_introns_coords):
     coords = map(lambda t: (t[1], t[0], t[2], t[3], t[4])
     if t[2] == '-' else (t[0], t[1], t[2], t[3], t[4]), exons_introns_coords)
     return coords

with open(infile) as f:
    dic_metadata = {}
    dic_genenames = {}  # {gene_id: gene_name}
    dic_gene_transcripts = {}  # {gene_id : {transcript_id: ([exons], {utr_type: [coords]})}}
    for line in f:
        if line.startswith('#'):  # metadata eg. #!genome-build
            line = line.strip().replace('#!', '').split()
            dic_metadata.setdefault(line[0], line[1])
        else:
            line = line.strip().split('\t')
            if line[2].lower().find('gene') != -1:
                line_from_gene_id = line[-1]  # gene_id "ENSXMAG00000022905"; ecc
                dic_from_parser = dic_gtf_column_val(line_from_gene_id)
                gene_id = dic_from_parser['gene_id']
                if dic_from_parser.get('gene_name') != None:
                    gene_name = dic_from_parser['gene_name']
                else:
                    gene_name = gene_id
                dic_genenames.setdefault(gene_id, gene_name)
                dic_gene_transcripts.setdefault(gene_id, ([], {}))  # nested dic initialization
                dic_gene_transcripts[gene_id][0]. \
                    append((int(line[3]), int(line[4]), line[6], line[0]))
            elif line[2].lower().find('exon') != -1:
                line_from_gene_id = line[-1]
                dic_from_parser = dic_gtf_column_val(line_from_gene_id)
                gene_id, transcript_id, exon_numb = dic_from_parser['gene_id'], \
                                                    dic_from_parser['transcript_id'], \
                                                    dic_from_parser['exon_number']

                dic_gene_transcripts[gene_id][1].setdefault(transcript_id, ([], {}))[0]. \
                    append((int(line[3]), int(line[4]), line[6], line[0], exon_numb))
            elif line[2].lower().find('utr') != -1:
                utr_type, line_from_gene_id = line[2], line[-1]
                dic_from_parser = dic_gtf_column_val(line_from_gene_id)
                gene_id, transcript_id, utr_type = dic_from_parser['gene_id'], \
                                                   dic_from_parser['transcript_id'], \
                                                   utr_type
                (dic_gene_transcripts[gene_id][1])[transcript_id][1].setdefault(utr_type, []). \
                    append((int(line[3]), int(line[4]), line[6], line[0]))

dic_gene_ids = {v: k for k, v in dic_genenames.items()}

header = ['Gene_name/Gene_ID', 'Gene_ID', 'Chrom', 'Coordinates', 'Transcript_id',
          'Exons', 'Introns','5_UTRs', '5_UTR_seqs', 'CAGEs', 'CAGEs_seq', '3_UTRs', '3_UTR_seqs']
header = '\t'.join(map(str, header))
print header
for key, value in dic_gene_transcripts.items():  # gene_id, splice_variants
    for k, v in value[1].items():  # transcript_coordinates [], utr {}
        if any(v[1]):  # if any utr {}
            gene_name_gene_id = dic_genenames.get(key)
            gene_id = dic_gene_ids.get(gene_name_gene_id)
            exons = v[0]
            introns = get_introns(exons)
            # exons_sequences = get_sequences(exons)
            # introns_sequences = get_sequences(introns)
            results_row = [gene_name_gene_id, gene_id, value[0][0][3], value[0][0],
                                   k, swap_coords(exons),
                                   swap_coords(introns), '', '', '', '', '', '']
            if v[1].get('five_prime_utr'):  # if 5'utr_
                five_prime_utr = v[1].get('five_prime_utr')
                five_prime_utr_seq = get_sequences(five_prime_utr)
                results_row[7] = five_prime_utr
                results_row[8] = five_prime_utr_seq
            if v[1].get('three_prime_utr'): #if 3' UTR
                three_prime_utr = v[1].get('three_prime_utr')
                three_prime_utr_seq = get_sequences(three_prime_utr)
                results_row[11] = three_prime_utr
                results_row[12] = three_prime_utr_seq
            print '\t'.join(map(str, results_row))
