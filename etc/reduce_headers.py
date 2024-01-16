#!/usr/bin/env python3
from Bio import SeqIO
import os
import gzip
import glob
import concurrent.futures

def reduce_headers(file_name):
    input_file = f"/home/glbrc.org/millican/repos/metagenome_snakemake/data/reads/{file_name}.fastq.gz"
    output_file = f"/home/glbrc.org/millican/repos/metagenome_snakemake/data/cor_reads/{file_name}.fastq.gz"
    with gzip.open(input_file, 'rt') as f, gzip.open(output_file, 'wt') as out:
        base_name = os.path.basename(f.name).replace(".fastq.gz", "")
        for record in SeqIO.parse(f, "fastq"):
            record.id = f"{base_name}_{record.id.split(':')[-1]}"
            record.description = f"{record.id} {record.description.split(' ')[-1]}"
            SeqIO.write(record, out, "fastq")    

def main():
    file_names = [os.path.basename(file.replace(".fastq.gz","")) for file in glob.glob("/home/glbrc.org/millican/repos/metagenome_snakemake/data/reads/*.fastq.gz")]
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(reduce_headers, file_names)
    
if __name__ == "__main__":
    main()