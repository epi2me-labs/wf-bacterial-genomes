#!/usr/bin/env extflow

nextflow.enable.dsl = 2



def helpMessage(){
    log.info """
Haploid Variant Analysis Workflow

Usage:
    nextflow run epi2me-labs/wf-hap-snp [options]

Options:
    --fastq             DIR     Path to FASTQ directory (required)
    --reference         FILE    Reference sequence FASTA file (required)
    --out_dir           DIR     Path for output (default: $params.out_dir)
    --medaka_model      STR     Medaka model name (default: $params.medaka_model)
    --run_prokka        BOOL    Run prokka on consensus sequence (default: $params.run_prokka)
    --prokka_opts       STR     Command-line arguments for prokka (default: $params.prokka_opts)

"""
}


process concatFastq {
    // concatenate fastq and fastq.gz in a dir

    label "containerCPU"
    cpus 1
    input:
        file "input"
    output:
        file "reads.fastq.gz"

    shell:
    '''
#!/usr/bin/env python
from glob import glob
import gzip
import itertools
import os
import pysam

# we use pysam just because it will read both fastq and fastq.gz
# and we don't have to worry about having a combination or not
with gzip.open("reads.fastq.gz", "wt") as fh:
    files = itertools.chain(
        glob("input/*.fastq"), glob("input/*.fastq.gz"))
    records = itertools.chain.from_iterable(
        pysam.FastxFile(fn) for fn in files) 
    for rec in records:
        annot = " {}".format(rec.comment) if rec.comment else ""
        qual = rec.quality if rec.quality else "+"*len(rec.sequence)
        fh.write("@{}{}\\n{}\\n+\\n{}\\n".format(rec.name, annot, rec.sequence, qual))
    '''
}


process overlapReads {
    // find overlaps of reads to reference

    label "containerCPU"
    cpus params.threads
    input:
        file reads
        file reference
    output:
        file "reads2ref.sam.gz"
        //tuple file("reads2ref.bam"), file("reads2ref.bam.csi")

    """
    # racon doesn't like bam 
    minimap2 -x map-ont --MD -a -t $task.cpus $reference $reads \
        | samtools sort -@ $task.cpus \
        | samtools view -h | bgzip > reads2ref.sam.gz
    """
}

process readStats {
    label "containerCPU"
    cpus 1
    input:
        file alignments
    output:
        file "readstats.txt"
        file "mapsummary.txt"
    """
    stats_from_bam -s mapsummary.txt -o readstats.txt $alignments
    """
}

process coverStats {
    label "containerCPU"
    cpus 2
    input:
        file alignments
    output:
        file "depth.txt"
    """
    # we need to convert sam to bam
    samtools view -@ $task.cpus -b --write-index $alignments -o reads.bam
    coverage_from_bam --one_file depth.txt --stride 100 reads.bam
    """
}

process scuffReference {
    // run racon from overlaps to obtain a consensus sequence

    label "containerCPU"
    cpus params.threads
    input:
        file reads
        file overlaps
        file reference
    output:
        file "racon.fa.gz"

    """
    racon --include-unpolished --no-trimming -q -1 -t $task.cpus \
        $reads $overlaps $reference | bgzip -c > racon.fa.gz
    """
}


process alignReadsToScuff {
    // align reads back to racon consensus sequence

    label "containerCPU"
    cpus params.threads
    input:
        file reads
        file ref
    output:
        tuple file("reads2scuffed.bam"), file("reads2scuffed.bam.bai")

    """
    minimap2 $ref $reads -x map-ont -t $task.cpus -a --secondary=no --MD -L | samtools sort --output-fmt BAM -o reads2scuffed.bam -@ $task.cpus
    samtools index reads2scuffed.bam -@ $task.cpus
    """
}


process splitRegions {
    // split the bam reference sequences into overlapping sub-regions

    label "containerCPU"
    cpus 1
    input:
        tuple file(bam), file(bai)
    output:
        stdout

    """
    #!/usr/bin/env python

    import itertools
    import medaka.common

    regions = itertools.chain.from_iterable(
        x.split($params.chunk_size, overlap=1000, fixed_size=False)
        for x in medaka.common.get_bam_regions("$bam"))
    for reg in regions:
        print(reg)
    """
}


// TODO: in a single GPU environment it would be better just
//       to use a single process for the whole bam file. Need
//       to read up on conditional channels

process medakaNetwork {
    // run medaka consensus for each region

    label "containerCPU"
    cpus 2
    input:
        tuple file(bam), file(bai)
        each reg
    output:
        file "consensus_probs.hdf"

    """
    medaka consensus $bam consensus_probs.hdf --region $reg
    """
}


process medakaConsensus {
    // gather all hdfs to create consolidated consensus sequence

    label "containerCPU"
    cpus params.threads
    input:
        file 'consensus_probs_*.hdf'
        file scuffed
    output:
        file "medaka_consensus.fasta"

    """
    medaka stitch --threads $task.cpus consensus_probs_*.hdf $scuffed medaka_consensus.fasta
    """
}


process medakaVCF {
    // create a VCF file by comparing the consensus to the reference

    label "containerCPU"
    cpus 1
    input:
        file fasta
        file reference
    output:
        file "medaka_consensus.vcf" 

    """
    medaka tools consensus2vcf $fasta $reference --out_prefix medaka_consensus
    """
}


process runProkka {
    // run prokka in a basic way on the consensus sequence

    label "prokka"
    cpus 1
    input:
        file "consensus.fasta"
    output:
        file "prokka_results"
    """
    prokka ${params.prokka_opts} --outdir prokka_results --prefix prokka consensus.fasta
    """
}


process makeReport {
    label "containerCPU"
    cpus 1
    input:
        file "depth.txt"
        file "read_summary.txt"
        file "align_summary.txt"
        file "variants.vcf"
    output:
        file "summary_report.html"
    """
    bcftools stats variants.vcf > variants.stats
    report.py depth.txt read_summary.txt align_summary.txt variants.stats summary_report.html
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


// modular workflow
workflow calling_pipeline {
    take:
        reads
        reference
    main:
        reads = concatFastq(reads)
        alignments = overlapReads(reads, reference)
        read_stats = readStats(alignments)
        depth_stats = coverStats(alignments)
        racon = scuffReference(reads, alignments, reference)
        aligns = alignReadsToScuff(reads, racon)
        regions = splitRegions(aligns).splitText()
        hdfs = medakaNetwork(aligns, regions)
        consensus = medakaConsensus(hdfs.collect(), racon)
        vcf = medakaVCF(consensus, reference)
        if (params.run_prokka) {
            prokka = runProkka(consensus)
        } else {
            prokka = Channel.empty()
        } 
        report = makeReport(depth_stats, read_stats[0], read_stats[1], vcf)
    emit:
        consensus
        vcf
        report
        prokka
}

// entrypoint workflow
workflow {
    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq or !params.reference) {
        helpMessage()
        println("")
        println("`--fastq` and `--reference` are required")
        exit 1
    }
    reads = channel.fromPath(params.fastq, type:'dir', checkIfExists:true)
    reference = channel.fromPath(params.reference) 
    results = calling_pipeline(reads, reference)
    output(results.consensus.concat(results.vcf, results.report, results.prokka))
}
