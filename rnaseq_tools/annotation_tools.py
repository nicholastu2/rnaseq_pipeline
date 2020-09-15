import re
from Bio import SeqIO
import pandas as pd
import os

def parseGtf_new(gtf_file, genome_fasta):
    """
        new gtf parser, written for NOIseq. eventually just use this for all uses
    """
    genome_dict = {}
    for seq_record in SeqIO.parse(genome_fasta, 'fasta'):
        genome_dict.setdefault(seq_record.name, seq_record.seq)

    print("...parsing genome")

    # gtfgenerator function
    with open(gtf_file, 'r') as gtf_file:
        parsed_gtf_dict = {}
        line_number = 1
        line = gtf_file.readline().strip().split("\t")
        print('...parsing gtf')
        while line:
            try:
                # extract features from gtf line
                chromosome = line[0]
                feature_type = line[2]
                start_coordinate = int(line[3])
                stop_coordinate = int(line[4])+1
                feature_length = stop_coordinate - start_coordinate
                strand = line[6]
                addtl_info_line = line[8].split(';')
                # extract gene id
                btwn_paren_regex = r"\"(.*)\""
                try:
                    gene_id = re.findall(btwn_paren_regex, addtl_info_line[0])[0]
                    if not gene_id.startswith('CNAG'):
                        raise ValueError('GeneIDNotCorrectFormat')
                except ValueError:
                    print('the extracted gene id did not start with CNAG on line %i of the gtf' %line_number)
                #extract biotype (the 'product' in the CNAG gtf this will be something like "hypothetical protein")
                try:
                    biotype = re.findall(btwn_paren_regex, addtl_info_line[3])[0]
                    if len(biotype) == 0:
                        raise ValueError('NoBiotype')
                except ValueError:
                    print('There is no biotype to extract on line %s of the gtf' %line_number)
                except IndexError:
                    if chromosome == 'NAT' or chromosome == 'G418':
                        biotype = "drug marker"
                    if chromosome == 'chrM':
                        biotype = 'hypothetical protein' # TODO: THIS I NEEDS TO BE FIXED IN THE GTF
                # add to dictionary
                gene_dict = {'chromosome': chromosome, 'strand': strand,'biotype': biotype, 'gc_content': 0}
                parsed_gtf_dict.setdefault(gene_id, gene_dict)
                # create entry in gene_id dictionary for cds_length
                if feature_type == 'CDS':
                    try:
                        parsed_gtf_dict[gene_id]['cds_length'] += feature_length
                    except KeyError:
                        # add 6 to the feature length for the stop and start codon (not included in CDS)
                        parsed_gtf_dict[gene_id].setdefault('cds_length', feature_length + 3)
                    finally:
                        cds_feature_sequence = genome_dict[chromosome][start_coordinate-1:stop_coordinate]
                        g_count = cds_feature_sequence.count("G")
                        g_count += cds_feature_sequence.count("g")
                        c_count = cds_feature_sequence.count("C")
                        c_count += cds_feature_sequence.count("c")
                        gc_count = g_count + c_count
                        parsed_gtf_dict[gene_id]['gc_content'] += gc_count
                if feature_type == 'start_codon':
                    parsed_gtf_dict[gene_id].setdefault('gene_start', start_coordinate)
                elif feature_type == 'stop_codon':
                    parsed_gtf_dict[gene_id].setdefault('gene_stop', stop_coordinate+2)

                line = gtf_file.readline().strip().split("\t")
                line_number +=1
            except IndexError:
                break
    return parsed_gtf_dict


def noiseqDataframe(parsed_gtf_dict, output_dir):
    """
        gtf data for R noiseq package. see https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
        for explanation of dataframes
    """
    mylength_dict = {}
    mygc_dict = {}
    mybiotypes_dict = {}
    mychroms_dict = {}
    for gene, gene_dict in parsed_gtf_dict.items():
        mylength_dict.setdefault(gene, gene_dict['cds_length'])
        mygc_dict.setdefault(gene, gene_dict['gc_content'] / float(gene_dict['cds_length']))
        mybiotypes_dict.setdefault(gene, gene_dict['biotype'])
        mychroms_dict.setdefault(gene_dict['chromosome'], [gene_dict['gene_start'], gene_dict['gene_stop']])

    mylength_df = pd.DataFrame.from_dict(mylength_dict, orient='index', columns=['length'])
    mygc_df = pd.DataFrame.from_dict(mygc_dict, orient='index', columns=['gc_content'])
    mybiotypes_df = pd.DataFrame.from_dict(mybiotypes_dict, orient='index', columns=['biotype'])
    mychroms_df = pd.DataFrame.from_dict(mychroms_dict, orient='index', columns=['geneStart', 'geneStop'])

    mylength_output_path = os.path.join(output_dir, 'mylength'+'.csv')
    mylength_df.to_csv(mylength_output_path)

    mygc_output_path = os.path.join(output_dir, 'mygc' + '.csv')
    mygc_df.to_csv(mygc_output_path)

    mybiotypes_output_path = os.path.join(output_dir, 'mybiotypes' + '.csv')
    mybiotypes_df.to_csv(mybiotypes_output_path)

    mychroms_output_path = os.path.join(output_dir, 'mychroms' + '.csv')
    mychroms_df.to_csv(mychroms_output_path)

# x = parseGtf_new('/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline/genome_files/KN99/crNeoKN99.gtf', '/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline/genome_files/KN99/crNeoKN99.fasta')
#
# noiseqDataframe(x, '/home/chase/code/cmatkhan/brentlab/scripts/old_crypto_eda_draft/data')

def parseGtf(gtf_file): # TODO: comment and docstring
    """
    Convert gtf into dictionary for filling bed fields
    :param gtf_file:
    :returns:
    """
    bed_dict = {}
    reader = open(gtf_file, 'r')
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


def parseGff3(gff_file): # TODO: comment and docstring
    """
    Convert gff3 into dictionary for filling bed fields
    :param gff_file:
    :returns:
    """
    bed_dict = {}
    reader = open(gff_file, 'r')
    for line in reader.readlines():
        if not line.startswith('#'):
            ## get info from each line
            line_split = line.split("\t")
            chrm, ltype, strand, annot = line_split[0], line_split[2], line_split[6], line_split[8]
            coords = [int(x) for x in line_split[3:5]]
            if ltype in 'gene':
                if annot.startswith('gene='): # TODO: this is a very ugly way of handling the fact that some lines are ID= and some are gene= in KN99 annote (marker is gene=)
                    gene_id = re.findall(r'gene=(.+?);', annot)[0]
                else:
                    gene_id = re.findall(r'ID=(.+?);', annot)[0]
                ## fill dictionary
                bed_dict[gene_id] = {'chrm': chrm, 'strand': strand, 'coords': [min(coords), max(coords)]}
    reader.close()
    return bed_dict