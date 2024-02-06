def calc_chunks(fasta):
    import gzip
    import math
    with gzip.open(fasta, "rt") as f:
        lines = f.readlines()
    return math.ceil((len(lines)/2)/12)

def get_mags(dir_path):
    import os
    files = os.listdir(dir_path)
    bin_list = [file.replace(".fna.gz", "") for file in files if file.replace(".fna.gz", "")]
    return bin_list

def create_tmpdir():
    import random
    import os
    import pickle
    with open("/home/glbrc.org/millican/repos/metagenome_snakemake/etc/adj-aml.pkl", 'rb') as f:
        adj, aml = pickle.load(f)
    temp_dir_base = "/home/glbrc.org/millican/TMPDIR"    # Replace with the base path for temporary directories
    # Construct the temporary directory path
    tmpdir = os.path.join(temp_dir_base, f"{random.choice(adj)}-{random.choice(aml)}")
    # Check if the directory exists, and find a new combination if it does
    while os.path.exists(tmpdir):
        tmpdir = os.path.join(temp_dir_base, f"{random.choice(adj)}-{random.choice(aml)}")
    # Once we find a combination that does not already exist
    # Create the temporary directory
    os.makedirs(tmpdir, exist_ok=True)
    return tmpdir

bin_names = get_mags(config["mag_directory"])

rule index:
    input:
        config["co_assembly"]
    output:
        "{data_dir}/index/co_assembly_index.1.bt2",
        "{data_dir}/index/co_assembly_index.2.bt2",
        "{data_dir}/index/co_assembly_index.3.bt2",
        "{data_dir}/index/co_assembly_index.4.bt2",
        "{data_dir}/index/co_assembly_index.rev.1.bt2",
        "{data_dir}/index/co_assembly_index.rev.2.bt2"
    params:
        "{data_dir}/index/co_assembly_index"
    threads: 12
    resources:
        mem_mb = 300000,
        tmpdir = create_tmpdir()
    conda:
        "read-mapping"
    shell:
        """
        bowtie2-build -f --threads {threads} {input} {params}
        """

rule map_reads:
    input:
        reads = "{data_dir}/reads/{sample_name}.fastq.gz",
        in1 = "{data_dir}/index/co_assembly_index.1.bt2",
        in2 = "{data_dir}/index/co_assembly_index.2.bt2",
        in3 = "{data_dir}/index/co_assembly_index.3.bt2",
        in4 = "{data_dir}/index/co_assembly_index.4.bt2",
        in5 = "{data_dir}/index/co_assembly_index.rev.1.bt2",
        in6 = "{data_dir}/index/co_assembly_index.rev.2.bt2"
    output:
        "{data_dir}/map/{sample_name}.sorted.bam"
    params:
        "{data_dir}/index/co_assembly_index"
    threads: 16
    resources:
        mem_mb = 300000,
        tmpdir = create_tmpdir()
    conda:
        "read-mapping"
    shell:
        """
        bowtie2 -x {params} --interleaved {input.reads} -p {threads} | samtools view -bS - | samtools sort -@ {threads} - > {output}
        """

rule index_bam:
    input:
        "{data_dir}/map/{sample_name}.sorted.bam"
    output:
        "{data_dir}/map/{sample_name}.sorted.bam.bai"
    threads: 16
    resources:
        mem_mb = 60000,
        tmpdir = create_tmpdir()
    conda:
        "read-mapping"
    shell:
        """
        samtools index -@ {threads} -o {output} {input}
        """

rule count_features:
    input:
        bam = "{data_dir}/map/{sample_name}.sorted.bam",
        idx = "{data_dir}/map/{sample_name}.sorted.bam.bai"
    output:
        "{data_dir}/count/{sample_name}.count"
    params:
        annotation = "{data_dir}/annotation/{sample_name}_functional_annotation.gff"
    threads: 16
    resources:
        mem_mb = 40000,
        tmpdir = create_tmpdir()
    conda:
        "read-mapping"
    shell:
        """
        featureCounts -T {threads} -p --countReadPairs --tmpDir {resources.tmpdir} -t CDS -g ID -F GFF -a {params.annotation} -o {output} {input.bam}
        """