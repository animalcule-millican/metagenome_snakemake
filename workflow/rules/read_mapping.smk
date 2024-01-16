configfile: "/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/config/read_assembly_config.yaml"
workdir: config["working_directory"]


def create_tmpdir():
    import random
    import os
    import pickle
    with open("/home/glbrc.org/millican/repos/metagenome_snakemake/etc/adj-aml.pkl", 'rb') as f:
        adj, aml = pickle.load(f)
    temp_dir_base = "/home/glbrc.org/millican/TMPDIR"
    tmpdir = os.path.join(temp_dir_base, f"{random.choice(adj)}-{random.choice(aml)}")
    while os.path.exists(tmpdir):
        tmpdir = os.path.join(temp_dir_base, f"{random.choice(adj)}-{random.choice(aml)}")
    os.makedirs(tmpdir, exist_ok=True)
    return tmpdir

rule index_coassembly:
    input:
        "{output_directory}/contigs/{project_name}.fna.gz"
    output:
        "{output_directory}/index/{project_name}_index.1.bt2",
        "{output_directory}/index/{project_name}_index.2.bt2",
        "{output_directory}/index/{project_name}_index.3.bt2",
        "{output_directory}/index/{project_name}_index.4.bt2",
        "{output_directory}/index/{project_name}_index.rev.1.bt2",
        "{output_directory}/index/{project_name}_index.rev.2.bt2"
    params:
        "{output_directory}/index/{project_name}_index"
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
        contig = "{output_directory}/contigs/{sample_name}.fna.gz",
        reads = "{output_directory}/reads/{sample_name}.ecc.fq.gz"
    output:
        "{output_directory}/bam/{sample_name}.sorted.bam"
    threads: 12
    resources:
        mem_mb = 10000,
        tmpdir = create_tmpdir()
    conda:
        "read-mapping"
    shell:
        """
        bowtie2-build -f --threads {threads} {input.contig} {resources.tmpdir}/{wildcards.sample_name}
        bowtie2 -x {resources.tmpdir}/{wildcards.sample_name} --interleaved {input.reads} -p {threads} | samtools view -bS - | samtools sort -@ {threads} - > {output}
        samtools index -@ {threads} {output}
        """

rule map_reads_coassembly:
    input:
        reads = "{output_directory}/reads/{sample_name}.ecc.fq.gz",
        in1 = "{output_directory}/index/{project_name}_index.1.bt2",
        in2 = "{output_directory}/index/{project_name}_index.2.bt2",
        in3 = "{output_directory}/index/{project_name}_index.3.bt2",
        in4 = "{output_directory}/index/{project_name}_index.4.bt2",
        in5 = "{output_directory}/index/{project_name}_index.rev.1.bt2",
        in6 = "{output_directory}/index/{project_name}_index.rev.2.bt2"
    output:
        "{data_dir}/map/{sample_name}.coassembly.sorted.bam"
    params:
        "{output_directory}/index/{project_name}_index"
    threads: 16
    resources:
        mem_mb = 300000,
        tmpdir = create_tmpdir()
    conda:
        "read-mapping"
    shell:
        """
        bowtie2 -x {params} --interleaved {input.reads} -p {threads} | samtools view -bS - | samtools sort -@ {threads} - > {output}
        samtools index -@ {threads} {output}
        rm -rf {resources.tmpdir}
        """



rule pileup_reads:
    shell:
        """
        pileup.sh in=aln.sam.gz out=cov.txt
        """

rule extract_unmapped:
    shell:
        """
        samtools view -u -f4 aln.sam.gz | samtools bam2fq -s unmapped.se.fq - > unmapped.pe.fq
        """