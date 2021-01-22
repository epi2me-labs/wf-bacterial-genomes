#!/usr/bin/env nextflow

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


// Both the input fastq and reference are used several times
Channel
    .fromPath(params.reads)
    .into { reads1; reads2; reads3 }

Channel
    .fromPath(params.reference)
    .into { reference1; reference2; reference3 }


process overlapReads {
    // find overlaps of reads to reference

    label "containerCPU"
    cpus params.threads

    input:
    file reads from reads1
    file reference from reference1

    output:
    file "reads2ref.paf" into reads2ref 

    """
    minimap2 -x map-ont -t $task.cpus $reference $reads > reads2ref.paf 
    """
}


process scuffReference {
    // run racon from overlaps to obtain a consensus sequence

    label "containerCPU"
    cpus params.threads

    input:
    file reads from reads2
    file paf from reads2ref
    file reference from reference2

    output:
    file "racon.fa.gz" into racon_consensus1, racon_consensus2

    """
    racon --include-unpolished --no-trimming -q -1 -t $task.cpus $reads $paf $reference | bgzip -c > racon.fa.gz
    """
}


process alignReadsToScuff {
    // align reads back to racon consensus sequence

    label "containerCPU"
    cpus params.threads

    input:
    file reads from reads3
    file ref from racon_consensus1
    
    output:
    file "reads2scuffed.bam" into medaka_bam1, medaka_bam2
    file "reads2scuffed.bam.bai" into medaka_bai1, medaka_bai2

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
    file bam from medaka_bam1
    file bai from medaka_bai1

    output:
    stdout into regions
    
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


regs = regions.splitText()


process medakaNetwork {
    // run medaka consensus for each region

    label "containerCPU"
    cpus 2

    input:
    file bam from medaka_bam2
    file bai from medaka_bai2
    each reg from regs

    output:
    file "consensus_probs.hdf" into medaka_hdf

    """
    medaka consensus $bam consensus_probs.hdf --region $reg
    """
}


// TODO: in a single GPU environment it would be better just
//       to use a single process for the whole bam file. Need
//       to read up on conditional channels


process medakaConsensus {
    // gather all hdfs to create consolidated consensus sequence

    label "containerCPU"
    cpus params.threads
    publishDir "${params.out_dir}", mode: 'copy', pattern: "medaka_consensus.fasta"

    input:
    file 'consensus_probs_*.hdf' from medaka_hdf.collect()
    file scuffed from racon_consensus2

    output:
    file "medaka_consensus.fasta" into medaka_fasta

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
    file fasta from medaka_fasta
    file reference from reference3

    output:
    file "medaka_consensus.vcf" 

    """
    medaka tools consensus2vcf $fasta $reference --out_prefix medaka_consensus
    """
}
