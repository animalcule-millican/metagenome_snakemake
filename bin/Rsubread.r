library("Rsubread")

# Parsed arguments
args <- Arg_Parser()
files = "/home/glbrc.org/millican/repos/metagenome_snakemake/data/map/Ga0501100.co_assembly.sorted.bam"
annot = "/home/glbrc.org/millican/projects/TMP_transfer_data/ptmp/adina/millican/oil/data/annotation/rhizo/rhizo-coassembly-annotation-feature-counting.gff"
tmpdir = "/home/glbrc.org/millican/TMPDIR"
threads = 4

# Count reads function
count_reads = function(files = files, annot = annot, tmpdir = tmpdir, threads = threads){
    require('Rsubread')
    df.fc = featureCounts(files,
        # annotation
        annot.inbuilt = NULL,
        annot.ext = annot,
        isGTFAnnotationFile = TRUE,
        GTF.featureType = "CDS",
        GTF.attrType = "ID",
        GTF.attrType.extra = c("product", "ko", "cog", "pfam", "tigrfam", "superfamily", "cath_funfam", "ec_number"),
        chrAliases = NULL,
        # level of summarization
        useMetaFeatures = TRUE,
        # overlap between reads and features
        allowMultiOverlap = FALSE,
        minOverlap = 1,
        fracOverlap = 0,
        fracOverlapFeature = 0,
        largestOverlap = FALSE,
        nonOverlap = NULL,
        nonOverlapFeature = NULL,
        # Read shift, extension and reduction
        readShiftType = "upstream",
        readShiftSize = 0,
        readExtension5 = 0,
        readExtension3 = 0,
        read2pos = NULL,
        # multi-mapping reads
        countMultiMappingReads = TRUE,
        # fractional counting
        fraction = FALSE,
        # long reads
        isLongRead = FALSE,
        # read filtering
        minMQS = 0,
        splitOnly = FALSE,
        nonSplitOnly = FALSE,
        primaryOnly = FALSE,
        ignoreDup = FALSE,
        # strandness
        strandSpecific = 0,
        # exon-exon junctions
        juncCounts = FALSE,
        genome = NULL,
        # parameters specific to paired end reads
        isPairedEnd = TRUE,
        countReadPairs = TRUE,
        requireBothEndsMapped = FALSE,
        checkFragLength = FALSE,
        minFragLength = 50,
        maxFragLength = 320,
        countChimericFragments = TRUE,
        autosort = TRUE,
        # number of CPU threads
        nthreads = threads,
        # read group
        byReadGroup = FALSE,
        # report assignment result for each read
        reportReads = NULL,
        reportReadsPath = NULL,
        # miscellaneous
        maxMOp = 10,
        tmpDir = tmpdir,
        verbose = FALSE)
    return(df.fc)
}

# Run Count Reads function
df.count = count_reads(files, annot, tmpdir, threads)

# Write output