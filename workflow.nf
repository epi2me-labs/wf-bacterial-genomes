#!/usr/bin/env nextflow

params.help = ""

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
    log.info ''

    return
}



reads_ch = Channel.fromPath(params.reads)
reference_ch = Channel.fromPath(params.reference)

process overlapReads {

    input:
    file reads from reads_ch
    file reference from reference_ch

    output:
    file "reads2ref.paf" into reads2ref 

    """
    echo $reference
    echo $reads
    minimap2 -x map-ont -t $task.cpus $reference $reads | bgzip -c > "reads2ref.paf" 
    """
}

process scuffReference {

    input:
    file reads from params.reads
    file paf from reads2ref
    file reference from params.reference

    output:
    file "racon.fa.gz" into racon_consensus

    """
    racon --include-unpolished --no-trimming -q -1 -t $task.cpus $reads $paf $reference | bgzip -c > racon.fa.gz
    """
}

process alignReadsToScuff {

    input:
    file reads from params.reads
    file ref from racon_consensus
    
    output:
    file "reads2scuffed.bam" into medaka_bam

    """
    minimap2 $ref $reads -x ont -t $task.cpus -a --secondary=no --MD -L | samtools sort --output-fmt BAM -o reads2scuffed.bam -@ $task.cpus
    samtools index reads2scuffed.bam -@ $task.cpus
    """
}

process medakaNetwork {

    input:
    file bam from medaka_bam

    output:
    file "consensus_probs.hdf" into medaka_hdf

    """
    medaka consensus $bam consensus_probs.hdf
    """
}

process medakaConsensus {

    input:
    file hdf from medaka_hdf

    output:
    file "medaka_consensus.fasta" into medaka_fasta

    """
    medaka stitch --fillgaps --threads $task.cpus $hdf medaka_consensus.fasta
    """
}

process medakaVCF {

    input:
    file fasta from medaka_fasta
    file reference from params.reference

    output:
    file "medaka.vcf" 
    stdout result

    """
    medaka tools consensus2vcf $fasta $reference --out_prefix medaka_consensus
    """
}

result.view { it.trim() }
