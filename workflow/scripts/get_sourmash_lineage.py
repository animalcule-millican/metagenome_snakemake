#!/usr/bin/env python3
import pickle
import pandas as pd
import sys
import taxopy
import concurrent.futures

def build_taxonomy_dict(dir_path):
    taxdb = taxopy.TaxDb(nodes_dmp=f"{dir_path}/nodes.dmp", names_dmp=f"{dir_path}/names.dmp", merged_dmp=f"{dir_path}/merged.dmp")
    return taxdb

def find_lineage(tax_list, taxdb):
    new_tax = {}
    acc = tax_list[0]
    taxid = tax_list[1]
    organism = tax_list[2]
    strain = tax_list[3]
    try:
        dat_dict = taxopy.Taxon(taxid, taxdb).rank_name_dictionary
        new_tax[acc] = {'ident': acc, 'taxid': taxid, 'superkingdom': dat_dict.get('superkingdom'), 'phylum': dat_dict.get('phylum'), 'class': dat_dict.get('class'), 'order': dat_dict.get('order'), 'family': dat_dict.get('family'), 'genus': dat_dict.get('genus'), 'species': organism, 'strain': strain}
    except taxopy.exceptions.TaxidError:
        return None
    return new_tax

def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    taxdmp = sys.argv[3]
    with open(input_file, 'rt') as f:
        data_dict = pickle.load(f)

    taxdb = build_taxonomy_dict(taxdmp)
    info_dict = {}
    lineage_dict = {}
    for key in data_dict.keys():
        info_dict[key] = {"accession": data_dict[key]['accession'], "taxid": data_dict[key]['taxid'], "organism": data_dict[key]['organism'], "strain": data_dict[key]['strain']}
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
        results = executor.map(find_lineage, info_dict.values(), taxdb)
    
    for result in results:
        if result is not None:
            lineage_dict.update(result)

    df = pd.DataFrame.from_dict(lineage_dict, orient='index')
    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    main()
            