configfile: "/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/config/read_assembly_config.yaml"
workdir: config["working_directory"]

def coassembly_input(input_directory):
    import glob
    import os
    file_list = glob.glob(f"{input_directory}/reads/*.norm.fq.gz")
    return ",".join(file_list)

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

def get_sample_names(input_directory):
    import glob
    import os
    sample_names = [os.path.basename(x).split(".")[0] for x in glob.glob(f"{input_directory}/*.fastq.gz")]
    return sample_names

sample_names = get_sample_names(config['sample_directory'])

wildcard_constraints:
    sample_name = "|".join(sample_names),
    project_name = config["project_name"]

rule all:
    input:
        expand("{output_directory}/contigs/{sample_name}.fna.gz", sample_name=config["sample_names"], output_directory=config["output_directory"]),
        expand("{output_directory}/contigs/{project_name}.fna.gz", project_name=config["project_name"], output_directory=config["output_directory"])

rule read_assembly:
    input:
        "{output_directory}/reads/{sample_name}.norm.fq.gz"
    output:
        "{output_directory}/contigs/{sample_name}.fna.gz"
    threads: 12
    resources:
        mem_mb = 300000
        tmpdir = create_tmpdir()
    conda:
        "assembly-tools"
    shell:
        """
        megahit --12 {input} -o {resources.tmpdir}/assembly --tmp-dir {resources.tmpdir} -t {threads} --presets meta-large
        gzip -c {resources.tmpdir}/assembly/final.contigs.fa > {output}
        """

rule coassembly:
    input:
        coassembly_input(config["output_directory"])
    output:
        "{output_directory}/contigs/{project_name}.fna.gz"
    threads: 8
    resources:
        mem_mb = 300000
        tmpdir = create_tmpdir()
    conda:
        "assembly-tools"
    shell:
        """
        megahit --12 {input} -o {resources.tmpdir}/assembly --tmp-dir {resources.tmpdir} -t {threads} --presets meta-large
        gzip -c {resources.tmpdir}/assembly/final.contigs.fa > {output}
        """

rule metabat2_depth:
    input:
        "{output_directory}/bam/{sample_name}.sorted.bam"
    output:
        "{output_directory}/depth/{sample_name}.sorted.bam.depth.txt"
    threads: 4
    resources:
        mem_mb = 30000,
        tmpdir = create_tmpdir()
    conda:
        "assembly-tools"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input}
        """
