#!/usr/bin/env python3
import taxoniq
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Convert NCBI genbank and refseq assembly summary to csv file with lineage, taxid and accesion information for sourmash classifications")
    parser.add_argument("-i", "--input", help="Input file from NCBI assembly summary")
    parser.add_argument("-o", "--output", help="Output file")
    args = parser.parse_args()
    return args

def get_lineage(acc, taxid, strain):
    lineage = {}
    tax_dict = {}
    tax_dict[taxid] = {}
    try:
        for t in taxoniq.Taxon(taxid).ranked_lineage:
            if t.rank.name in ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]:
                tax_dict[taxid].update({t.rank.name: t.scientific_name})
        lineage[acc] = {
            "ident": acc, 
            "taxid": taxid, 
            "superkingdom": tax_dict[taxid].get("superkingdom", ""), 
            "phylum": tax_dict[taxid].get("phylum", ""), 
            "class": tax_dict[taxid].get("class", ""), 
            "order": tax_dict[taxid].get("order", ""), 
            "family": tax_dict[taxid].get("family", ""), 
            "genus": tax_dict[taxid].get("genus", ""), 
            "species": tax_dict[taxid].get("species", ""), 
            "strain": strain
        }
    except KeyError:
        print(f"Taxid {taxid} not found in the database.")
    return lineage

def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    lineage_dict = {}
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("#") or line.startswith("assembly_accession"):
                continue
            line = line.strip().split("\t")
            acc = line[0]
            taxid = line[5]
            if line[8] == "na":
                strain = ""
            else:
                strain = f"{line[7]} {line[8].replace('strain=', '')}"
            lineage_dict.update(get_lineage(acc, taxid, strain))
    
    df = pd.DataFrame.from_dict(lineage_dict, orient="index")
    df.to_csv(output_file, sep=",", index=False)

if __name__ == "__main__":
    main()