configfile: "/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/modules/taxonomy/taxonomy_config.yml"
workdir: config["working_directory"]

wildcard_constraints:
    taxa = "|".join(config["taxa"]),
    database = "|".join(config["databases"])

rule all:
    input:


rule get_reports:
    input:
        taxa = "{taxa}",
        database = "{database}"
    output:
        "{reference_directory}/assembly_reports/{taxa}-{database}_assembly_summary.pkl"
    threads: 1
    resources:
        mem_mb = 4000
    conda:
        "prepare-reference"
    shell:
        """
        scripts/pickle_assembly_reports.py {input.taxa} {output} {input.database}
        """

rule get_ncbi_files:
    input:
        "{reference_directory}/pkl/{taxa}-{database}_assembly_summary.pkl"
    output:
        "{reference_directory}/pkl/ncbi_{taxa}-{database}_genome_info.pkl"
    params:
        ref = config["reference_directory"]
    threads: 6
    resources:
        mem_mb = 18000
    conda:
        "prepare-reference"
    shell:
        """
        scripts/get_genome_files.py {input} {params.ref} {output} {wildcards.taxa}
        """

rule parse_ncbi_taxonomy:
    input:
        "{reference_directory}/pkl/{taxa}-{database}_assembly_summary.pkl"
    output:
        "{reference_directory}/taxonomy/{taxa}-{database}_ncbi_taxonomy.csv"
    threads: 1
    resources:
        mem_mb = 8000
    conda:
        "prepare-reference"
    shell:
        """
        scripts/parse_ncbi_taxonomy_taxidTools.py {input} {output}
        """

rule join_taxonomy:
    input:
        expand("{reference_directory}/taxonomy/{{taxa}}-{database}_ncbi_taxonomy.csv", refdir = config["reference_directory"], taxa = config["taxa"], database = config["database"]),
    output:
        "{reference_directory}/taxonomy/{taxa}_reference_taxonomy.csv"
    params:
        refdir = config["reference_directory"]
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "prepare-reference"
    shell:
        """
        scripts/join-taxonomy.py {params.refdir} {input}
        """

rule taxonomy_db:
    input:
        expand("{reference_directory}/taxonomy/{{taxa}}-{database}_ncbi_taxonomy.csv", refdir = config["reference_directory"], database = config["database"]),
    output:
        "{reference_directory}/lineage/{taxa}_reference_taxonomy.db"
    threads: 1
    resources:
        mem_mb=24000
    conda:
        "branchwater"
    shell:
        """
        sourmash tax prepare --taxonomy {input} -o {output}
        """

rule sketch_references:
    input:
        expand("{reference_directory}/pkl/ncbi_{{taxa}}-{database}_genome_info.pkl", refdir = config["reference_directory"], database = config["database"])
    output:
        "{reference_directory}/sketch/{taxa}_sketch.zip"
    conda:
        "branchwater"
    threads: 12
    resources:
        mem_mb=30000
    shell:
        """
        sourmash scripts manysketch -p k=21,K=31,k=51,scaled=1000,abund -c {threads} -o {output} {input.genome_list}
        """

rule reference_sketch_index:
    input:
        sketchfile = "{reference_directory}/sketch/{taxa}_sketch.zip",
        kmer = "{kmer}"
    output:
        "{reference_directory}/sketch/index/{taxa}_sketch.k{kmer}.idx"
    threads: 12
    resources:
        mem_mb=60000
    conda:
        "branchwater"
    shell:
        """
        sourmash scripts index --ksize {wildcards.kmer} --scaled 1000 --output {output} --cores {threads} {input.sketchfile}
        """