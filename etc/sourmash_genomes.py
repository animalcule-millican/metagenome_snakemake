#!/usr/bin/env python3
import urllib.request
import pandas as pd
import sys
import taxopy
import os

def download_genome(url, genome_path):
    file_path = f"{genome_path}/{url.split('/')[-1]}"
    if not os.path.exists(file_path):
        urllib.request.urlretrieve(url, file_path)
        return file_path
    else:
        print(f"File {file_path} already exists, skipping download.")
    return None

def build_taxonomy_dict(dir_path = "/home/glbrc.org/millican/ref_db/taxdmp"):
    taxdb = taxopy.TaxDb(nodes_dmp=f"{dir_path}/nodes.dmp", names_dmp=f"{dir_path}/names.dmp", merged_dmp=f"{dir_path}/merged.dmp")
    return taxdb

def find_lineage(taxid, acc, taxdb, organism, strain):
    new_tax = {}
    taxid = int(taxid)
    try:
        dat_dict = taxopy.Taxon(taxid, taxdb).rank_name_dictionary
        new_tax[acc] = {'ident': acc, 'taxid': taxid, 'superkingdom': dat_dict.get('superkingdom'), 'phylum': dat_dict.get('phylum'), 'class': dat_dict.get('class'), 'order': dat_dict.get('order'), 'family': dat_dict.get('family'), 'genus': dat_dict.get('genus'), 'species': organism, 'strain': strain}
    except taxopy.exceptions.TaxidError as e:
        print(e)
        print(taxid)
        return new_tax
    return new_tax

def main():
    taxa = sys.argv[1]
    genome_path = f"/home/glbrc.org/millican/ref_db/reference_genomes/genomes/{taxa}"
    input_file = f"/home/glbrc.org/millican/ref_db/reference_genomes/assembly_summary_{taxa}.txt"

    if not os.path.exists(genome_path):
        os.makedirs(genome_path)
    lineage_dict = {}
    taxdb = build_taxonomy_dict(dir_path = "/home/glbrc.org/millican/ref_db/taxdmp")
    with open(input_file, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            row = line.strip().split('\t')
            accession = row[0]
            taxid = row[6]
            organism = row[7]
            if row[8] != 'na':
                strain = row[8].split('=')[1]
            elif row[8] == 'na':
                strain = ''
            ftp = f"{row[19]}/{row[19].split('/')[-1]}_genomic.fna.gz"
            dw_ftp = download_genome(ftp, genome_path)           
            out = find_lineage(taxid, accession, taxdb, organism, strain)
            lineage_dict.update(out)
            

    df = pd.DataFrame.from_dict(lineage_dict, orient='index')
    df.to_csv(f"/home/glbrc.org/millican/ref_db/reference_genomes/{taxa}_lineage.csv", index=False)

if __name__ == "__main__":
    main()
            