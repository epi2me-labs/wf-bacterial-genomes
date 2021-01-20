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


Channel
    .fromPath(params.reads)
    .into { reads1; reads2; reads3 }

Channel
    .fromPath(params.reference)
    .into { reference1; reference2; reference3 }


process overlapReads {

    label "containerCPU"
    cpus params.threads

    input:
    file reads from reads1
    file reference from reference1

    output:
    file "reads2ref.paf" into reads2ref 

    """
    echo $reference
    echo $reads
    echo "CPUS: "$task.cpus
    minimap2 -x map-ont -t $task.cpus $reference $reads > "reads2ref.paf" 
    """
}


process scuffReference {

    label "containerCPU"
    cpus params.threads

    input:
    file reads from reads2
    file paf from reads2ref
    file reference from reference2

    output:
    file "racon.fa.gz" into racon_consensus1, racon_consensus2

    """
    echo "CPUS: "$task.cpus
    racon --include-unpolished --no-trimming -q -1 -t $task.cpus $reads $paf $reference | bgzip -c > racon.fa.gz
    """
}

process alignReadsToScuff {

    label "containerCPU"
    cpus params.threads

    input:
    file reads from reads3
    file ref from racon_consensus1
    
    output:
    file "reads2scuffed.bam" into medaka_bam
    file "reads2scuffed.bam.bai" into medaka_bai

    """
    minimap2 $ref $reads -x map-ont -t $task.cpus -a --secondary=no --MD -L | samtools sort --output-fmt BAM -o reads2scuffed.bam -@ $task.cpus
    samtools index reads2scuffed.bam -@ $task.cpus
    """
}

process medakaNetwork {

    label "containerGPU"
    cpus 2

    input:
    file bam from medaka_bam
    file bai from medaka_bai

    output:
    file "consensus_probs.hdf" into medaka_hdf

    """
    medaka consensus $bam consensus_probs.hdf
    """
}

process medakaConsensus {

    label "containerCPU"
    cpus params.threads
    publishDir "${params.out_dir}", mode: 'copy', pattern: "medaka_consensus.fasta"

    input:
    file hdf from medaka_hdf
    file scuffed from racon_consensus2

    output:
    file "medaka_consensus.fasta" into medaka_fasta

    """
    medaka stitch --threads $task.cpus $hdf $scuffed medaka_consensus.fasta
    """
}

process medakaVCF {

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
