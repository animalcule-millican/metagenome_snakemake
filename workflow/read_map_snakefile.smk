
configfile: "/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/config.yml"
workdir: config["working_directory"]
def get_names():
    import glob
    import os
    names_list = [os.path.basename(x).replace(".fastq.gz", "") for x in glob.glob("/home/glbrc.org/millican/repos/metagenome_snakemake/data/reads/*.fastq.gz")]
    return names_list


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

names = get_names()

wildcard_constraints:
    sample_name = "|".join(names),
    data_dir = config["data_directory"]

rule all:
    input:
        expand("{data_dir}/count/{sample_name}.co_assembly.count", sample_name=names, data_dir=config["data_directory"])


rule map_reads_coassembly:
    input:
        reads = "{data_dir}/reads/{sample_name}.fastq.gz",
        in1 = "{data_dir}/index/co_assembly_index.1.bt2l",
        in2 = "{data_dir}/index/co_assembly_index.2.bt2l",
        in3 = "{data_dir}/index/co_assembly_index.3.bt2l",
        in4 = "{data_dir}/index/co_assembly_index.4.bt2l",
        in5 = "{data_dir}/index/co_assembly_index.rev.1.bt2l",
        in6 = "{data_dir}/index/co_assembly_index.rev.2.bt2l"
    output:
        "{data_dir}/map/{sample_name}.co_assembly.sorted.bam"
    params:
        idx = "{data_dir}/index/co_assembly_index",
        tmpdir = create_tmpdir()
    threads: 8
    resources:
        mem_mb = 500000
    conda:
        "read-mapping"
    shell:
        """
        export TMPDIR={params.tmpdir}
        bowtie2 -x {params.idx} --interleaved {input.reads} -p {threads} | samtools view -bS - | samtools sort -@ {threads} - > {output}
        rm -rf {params.tmpdir}
        """

rule index_maps:
    input:
        "{data_dir}/map/{sample_name}.co_assembly.sorted.bam"
    output:
        "{data_dir}/map/{sample_name}.co_assembly.sorted.bai"
    params:
        tmpdir = create_tmpdir()
    threads: 4
    resources:
        mem_mb = 100000
    conda:
        "read-mapping"
    shell:
        """
        export TMPDIR={params.tmpdir}
        samtools index -@ {threads} -o {output} {input}
        rm -rf {params.tmpdir}
        """

rule count_features:
    input:
        bam = "{data_dir}/map/{sample_name}.co_assembly.sorted.bam",
        bai = "{data_dir}/map/{sample_name}.co_assembly.sorted.bai"
    output:
        "{data_dir}/count/{sample_name}.co_assembly.count"
    params:
        gff = config["gff"],
        tmpdir = create_tmpdir()
    threads: 1
    resources:
        mem_mb = 200000
    conda:
        "read-mapping"
    shell:
        """
        #export TMPDIR={params.tmpdir}
        #scripts/count_reads.r -i {input.bam} -a {params.gff} -t {params.tmpdir} -n {threads} -o {output}
        export TMPDIR={params.tmpdir}
        htseq-count -s no -r pos -t CDS -i ID -f bam {input.bam} {params.gff} > {output}
        rm -rf {params.tmpdir}
        """