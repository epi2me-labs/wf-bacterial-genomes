#!/usr/bin/env extflow

nextflow.enable.dsl = 2


params.help = ""
params.threads = 1

if(params.help) {
    log.info ''
    log.info 'Haploid Variant Calling'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run workflow.nf [options]'
    log.info ''
    log.info 'Script Options: '
    log.info '    --reads        FILE    Path to FASTQ file'
    log.info '    --reference    FILE    Path to reference genome'
    log.info '    --out_dir      PATH    Output directory'
    log.info '    --threads      INT     Number of CPU threads to use'
    log.info ''

    return
}


process overlapReads {
    // find overlaps of reads to reference

    label "containerCPU"
    cpus params.threads

    input:
    file reads
    file reference

    output:
    file "reads2ref.paf" 

    """
    minimap2 -x map-ont -t $task.cpus $reference $reads > reads2ref.paf 
    """
}


process scuffReference {
    // run racon from overlaps to obtain a consensus sequence

    label "containerCPU"
    cpus params.threads

    input:
    file reads
    file paf
    file reference

    output:
    file "racon.fa.gz"

    """
    racon --include-unpolished --no-trimming -q -1 -t $task.cpus $reads $paf $reference | bgzip -c > racon.fa.gz
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
        x.split(int(1e6), overlap=1000, fixed_size=False)
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
    publishDir "${params.out_dir}", mode: 'copy', pattern: "medaka_consensus.fasta"

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
    publishDir "${params.out_dir}", mode: 'copy', pattern: "medaka_consensus.vcf"

    input:
    file fasta
    file reference

    output:
    file "medaka_consensus.vcf" 

    """
    medaka tools consensus2vcf $fasta $reference --out_prefix medaka_consensus
    """
}


workflow calling_pipeline {
    take:
        reads
        reference
    main:
        overlaps = overlapReads(reads, reference)
        racon = scuffReference(reads, overlaps, reference)
        aligns = alignReadsToScuff(reads, racon)
        regions = splitRegions(aligns).splitText()
        hdfs = medakaNetwork(aligns, regions)
        consensus = medakaConsensus(hdfs.collect(), racon)
        vcf = medakaVCF(consensus, reference)
    emit:
        consensus
        vcf
}


workflow {
    reads = channel.fromPath(params.reads)
    reference = channel.fromPath(params.reference) 
    calling_pipeline(reads, reference)
    // TODO: how to we publish files from here rather than the processes?
}
