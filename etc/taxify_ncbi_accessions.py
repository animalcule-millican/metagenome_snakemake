#!/usr/bin/env python3
import pandas as pd
import sys
import os
import pickle
import taxopy
import subprocess
import concurrent.futures



def download_genomes(lin_dict, key, out_dir, taxa):
    sour_script = "/home/glbrc.org/millican/repos/metagenome_snakemake/etc/sourmash-get-genomes.sh"
    acc = lin_dict[key]['ident']
    cmd = ['bash', sour_script, os.path.join(out_dir, taxa), acc]
    subprocess.run(cmd)

def build_taxonomy_dict(dir_path = "/home/glbrc.org/millican/ref_db/taxdmp"):
    taxdb = taxopy.TaxDb(nodes_dmp=f"{dir_path}/nodes.dmp", names_dmp=f"{dir_path}/names.dmp", merged_dmp=f"{dir_path}/merged.dmp")
    return taxdb

def find_lineage(taxid_dict, taxdb, tax_keys):
    new_tax = {}
    invalid_list = []
    for key in taxid_dict.keys():
        taxid = int(taxid_dict[key]['taxid'].strip())
        try:
            dat_dict = taxopy.Taxon(taxid, taxdb).rank_name_dictionary
            if dat_dict.get('strain') == None:
                strain = ''
            else:
                strain = dat_dict.get('strain')
            
            new_tax[key] = {'ident': key, 'taxid': taxid, 'superkingdom': dat_dict.get('superkingdom'), 'phylum': dat_dict.get('phylum'), 'class': dat_dict.get('class'), 'order': dat_dict.get('order'), 'family': dat_dict.get('family'), 'genus': dat_dict.get('genus'), 'species': dat_dict.get('species'), 'strain': strain}
            
        except taxopy.exceptions.TaxidError:
            invalid_list.append(taxid)
    return new_tax, invalid_list

def parse_assembly_summary(genbank_file, refseq_file):
    info_dict = {}
    acc_set = []
    with open(genbank_file, 'r') as f:
        for i, line in enumerate(f):
            if i < 2:
                continue
            row = line.strip().split('\t')
            info_dict[row[0]] = {"accession": row[0], 'taxid': row[5]}
            acc_set.append(row[0].replace("GCA", ""))
    acc_list = set(acc_set)
    with open(refseq_file, 'r') as f:
        for i, line in enumerate(f):
            if i < 2:
                continue
            row = line.strip().split('\t')
            if not row[0].replace("GCF", "") in acc_list:
                info_dict[row[0]] = {"accession": row[0], 'taxid': row[5]}
    return info_dict

def main():
    out_dir = "/home/glbrc.org/millican/ref_db/reference_genomes/genomes"
    taxa = sys.argv[1]
    taxdb = build_taxonomy_dict(dir_path = "/home/glbrc.org/millican/ref_db/taxdmp")
    tax_keys = ['superkingdom','phylum','class','order','family','genus','species','strain']
    invalid_taxids = []
    genbank_file = f"/home/glbrc.org/millican/ref_db/sourDB/util-files/assembly_summary_genbank_{taxa}.txt"
    refseq_file = f"/home/glbrc.org/millican/ref_db/sourDB/util-files/assembly_summary_refseq_{taxa}.txt"
    data_dict = parse_assembly_summary(genbank_file, refseq_file)
    lin_dict, invalid = find_lineage(data_dict, taxdb, tax_keys)
    with open(f"/home/glbrc.org/millican/ref_db/sourDB/util-files/{taxa}_dict.pkl", 'wb') as f:
        pickle.dump(lin_dict, f, protocol = pickle.HIGHEST_PROTOCOL)
    invalid_taxids.append(invalid)
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for key in lin_dict.keys():
            executor.submit(download_genomes, lin_dict, key, out_dir, taxa)
    
    df = pd.DataFrame.from_dict(lin_dict, orient = 'index')
    df.to_csv(f"/home/glbrc.org/millican/ref_db/reference_genomes/taxonomy/{taxa}_lineage.csv", index = False)
    genbank_file = ""
    refseq_file = ""
    data_dict = {}
    lin_dict = {}
    invalid = []
                
      
if __name__ == "__main__":
    main()