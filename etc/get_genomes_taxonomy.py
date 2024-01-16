#!/usr/bin/env python3
import urllib.request
import pandas as pd
import sys
import taxopy


def download_genome(url, genome_path):
    file_path = f"{genome_path}/{url.split('/')[-1]}"
    urllib.request.urlretrieve(url, file_path)

def build_taxonomy_dict(dir_path = "/home/glbrc.org/millican/ref_db/taxdmp"):
    taxdb = taxopy.TaxDb(nodes_dmp=f"{dir_path}/nodes.dmp", names_dmp=f"{dir_path}/names.dmp", merged_dmp=f"{dir_path}/merged.dmp")
    return taxdb

def find_lineage(taxid, acc, taxdb):
    new_tax = {}
    try:
        dat_dict = taxopy.Taxon(taxid, taxdb).rank_name_dictionary
        if dat_dict.get('strain') == None:
            strain = ''
        else:
            strain = dat_dict.get('strain')
        new_tax[acc] = {'ident': acc, 'taxid': taxid, 'superkingdom': dat_dict.get('superkingdom'), 'phylum': dat_dict.get('phylum'), 'class': dat_dict.get('class'), 'order': dat_dict.get('order'), 'family': dat_dict.get('family'), 'genus': dat_dict.get('genus'), 'species': dat_dict.get('species'), 'strain': strain}
    except taxopy.exceptions.TaxidError:
        return None
    return new_tax

def get_taxid_acc(url):
    with urllib.request.urlopen(url) as response:
        content = response.read().decode("utf-8").splitlines()
        acc = None
        taxid = None
        for line in content:
            if "Taxid:" in line:
                line = line.strip()
                pos = line.find("Taxid:")
                assert pos >= 0
                pos += len("Taxid:")
                taxid = line[pos:]
                taxid = taxid.strip()
            if "assembly accession:" in line:
                line = line.strip()
                pos = line.find("assembly accession:")
                assert pos >= 0
                pos += len("assembly accession:")
                acc = line[pos:]
                acc = acc.strip()
            if taxid and acc is not None:
                return taxid, acc

def get_ftp_urls(input_file):
    info_dict = {}
    genome_list = []
    report_list = []
    with open(input_file, 'r') as f:
        for i, line in enumerate(f):
            if i < 2:
                continue
            row = line.strip().split('\t')
            info_dict[row[0]] = {"accession": row[0], 'taxid': row[5], 'assembly_name': row[15], 'ftp_path': row[19]}
    for key in info_dict.keys():
        rootname = info_dict[key]['ftp_path'].split('/')[-1]
        report_url = f"{info_dict[key]['ftp_path']}/{rootname}_assembly_report.txt"
        genome_url = f"{info_dict[key]['ftp_path']}/{rootname}_genomic.fna.gz"
        genome_list.append(genome_url)
        report_list.append(report_url)
    return genome_list, report_list

def merge_summary_files(genbank_file, refseq_file, taxa):
    # Read the first TSV file
    df1 = pd.read_csv(genbank_file, sep='\t')
    # Read the second TSV file
    df2 = pd.read_csv(refseq_file, sep='\t')
    # Merge the two dataframes based on common headers
    common_columns = df1.columns.intersection(df2.columns)
    # Merge on common columns
    merged_df = pd.merge(df1, df2, on=common_columns)
    # Write the merged dataframe to a new TSV file
    outfile = f"/home/glbrc.org/millican/ref_db/sourDB/util-files/assembly_summary_{taxa}.txt"
    merged_df.to_csv(outfile, sep='\t', index=False)

def main():
    taxa = sys.argv[1]
    taxdb = build_taxonomy_dict(dir_path = "/home/glbrc.org/millican/ref_db/taxdmp")
    tax_dict = {}
    genome_path = f"/home/glbrc.org/millican/ref_db/reference_genomes/genomes/{taxa}"
    genbank_file = f"/home/glbrc.org/millican/ref_db/sourDB/util-files/assembly_summary_genbank_{taxa}.txt"
    refseq_file = f"/home/glbrc.org/millican/ref_db/sourDB/util-files/assembly_summary_refseq_{taxa}.txt"
    merge_summary_files(genbank_file, refseq_file, taxa)
    summary_file = f"/home/glbrc.org/millican/ref_db/sourDB/util-files/assembly_summary_{taxa}.txt"
    genome_list, report_list = get_ftp_urls(summary_file)
    
    for url in report_list:
        taxid, acc = get_taxid_acc(url)
        tax_dict.update(find_lineage(taxid, acc, taxdb))
    
    for url in genome_list:
        download_genome(url, genome_path)
    
    df = pd.DataFrame.from_dict(tax_dict, orient = 'index')
    df.to_csv(f"/home/glbrc.org/millican/ref_db/reference_genomes/taxonomy/{taxa}_lineage.csv", index = False)
    
if __name__ == "__main__":
    main()