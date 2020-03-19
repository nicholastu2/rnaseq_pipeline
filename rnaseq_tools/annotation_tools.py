
def parse_gtf(filename):
    """
    Convert gtf into dictionary for filling bed fields
    """
    bed_dict = {}
    reader = open(filename, 'r')
    for line in reader.readlines():
        if not line.startswith('#'):
            ## get info from each line
            line_split = line.split("\t")
            chrm, ltype, strand, annot = line_split[0], line_split[2], line_split[6], line_split[8]
            coords = [int(x) for x in line_split[3:5]]
            gene_id = re.findall(r'gene_id "(.+?)";', annot)[0]
            ## fill dictionary
            if gene_id not in bed_dict.keys():
                bed_dict[gene_id] = {'chrm': chrm, 'strand': strand, 'coords': [0, 0]}
            if ltype in 'start_codon':
                if strand == '+':
                    bed_dict[gene_id]['coords'][0] = min(coords)
                elif strand == '-':
                    bed_dict[gene_id]['coords'][1] = max(coords)
            elif ltype in 'stop_codon':
                if strand == '+':
                    bed_dict[gene_id]['coords'][1] = max(coords)
                elif strand == '-':
                    bed_dict[gene_id]['coords'][0] = min(coords)
    reader.close()
    return bed_dict


def parse_gff3(filename):
    """
    Convert gff3 into dictionary for filling bed fields
    """
    bed_dict = {}
    reader = open(filename, 'r')
    for line in reader.readlines():
        if not line.startswith('#'):
            ## get info from each line
            line_split = line.split("\t")
            chrm, ltype, strand, annot = line_split[0], line_split[2], line_split[6], line_split[8]
            coords = [int(x) for x in line_split[3:5]]
            if ltype in 'gene':
                gene_id = re.findall(r'Name=(.+?);', annot)[0]
                ## fill dictionary
                bed_dict[gene_id] = {'chrm': chrm, 'strand': strand, 'coords': [min(coords), max(coords)]}
    reader.close()
    return bed_dict