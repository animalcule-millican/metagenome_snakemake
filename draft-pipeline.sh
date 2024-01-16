#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate read-mapping
source random_directory.sh
oldPERL5LIB=/opt/bifxapps/perl516/share/perl5:/opt/bifxapps/perl516/lib64/perl5:/opt/bifxapps/perl516/lib/perl5
export PERL5LIB=/home/glbrc.org/millican/mambaforge/envs/read-mapping/lib/perl5
export samtools=/home/glbrc.org/millican/mambaforge/envs/read-mapping/bin/samtools
export bowtie=/home/glbrc.org/millican/mambaforge/envs/read-mapping/bin/bowtie2
export index=/home/glbrc.org/millican/projects/oil/data/map-index/rhizo/rhizo-coassembly
bowtie2-build  --threads $3 -f $1 $2

$bowtie --threads $3 --mm -x $index --interleaved $1 > $2

rm -r $TMPDIR
export PERL5LIB=$oldPERL5LIB

/home/glbrc.org/millican/projects/TMP_transfer_data/ptmp/adina/millican/oil/data/assembly/rhizo/rhizo-coassembly.contigs.fna.gz
##### Rhizosphere Co-Assembly 
# annotation /home/glbrc.org/millican/projects/TMP_transfer_data/ptmp/adina/millican/oil/data/annotation/rhizo/Ga0531775_functional_annotation.gff
# faa /home/glbrc.org/millican/projects/TMP_transfer_data/ptmp/adina/millican/oil/data/annotation/rhizo/rhizo-coassembly.faa
salmon index -t metagG.ffn -i transcript_index --type quasi -k 31
salmon quant -i transcript_index -1 $BASE$tail1 -2 $BASE$tail2 -o $BASE.quant --meta -g $GFF
salmon quant --alignments $1 --targets --output $2 --threads $cpu --mappingCacheMemoryLimit 100000000
#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate read-mapping
source random_directory.sh
export samtools=/home/glbrc.org/millican/mambaforge/envs/read-mapping/bin/samtools
$samtools view -bS $1 | $samtools sort -@ $3 - > $2

#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate read-mapping
source random_directory.sh
export samtools=/home/glbrc.org/millican/mambaforge/envs/read-mapping/bin/samtools
$samtools index -@ 4 $OUT

#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate read-mapping
source random_directory.sh
export samtools=/home/glbrc.org/millican/mambaforge/envs/read-mapping/bin/samtools
$samtools coverage $OUT > $DIR/data/stats/samples/mags/${SAMPLE_NAME}-${1}.coverage.txt

#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate read-mapping
source random_directory.sh
export samtools=/home/glbrc.org/millican/mambaforge/envs/read-mapping/bin/samtools
$samtools depth -@ 4 $OUT > $DIR/data/stats/samples/mags/${SAMPLE_NAME}-${1}.depth.txt

#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate read-mapping
source random_directory.sh
export samtools=/home/glbrc.org/millican/mambaforge/envs/read-mapping/bin/samtools
$samtools idxstats -@ 4 $OUT > $DIR/data/stats/samples/mags/${SAMPLE_NAME}-${1}.idxstats.txt

#!/bin/bash
source /mnt/bigdata/linuxhome/millican/.bashrc
conda activate thats-so-meta
featureCounts -T 2 -p --countReadPairs -t CDS -g "${feature}" -f -a "$GFF" -o "$COUNTS/${NAME}_${feature}.counts.txt" "$INPUT"



### Taxonomy assignment for reads