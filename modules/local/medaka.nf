process medakaHdf {
    // run medaka consensus for an individual region
    label "medaka"
    cpus 2
    // medaka rarely uses more than 8 GB, but sometimes it does happen
    memory { task.attempt == 1 ? "8 GB" : "15 GB" }
    errorStrategy { task.exitStatus == 137 ? "retry" : "terminate" }
    maxRetries 1
    input:
        tuple val(meta),
            path("align.bam"),
            path("align.bam.bai"),
            val(region),
            val(basecall_model)
        val type
    output:
        tuple val(meta), path("*consensus_probs.hdf")
    script:
    assert type in ["consensus", "variant"]
    """
    medaka --version
    echo ${basecall_model}
    medaka consensus align.bam "${meta.alias}.consensus_probs.hdf" \
        --threads 2 --regions "${region}" --model ${basecall_model}:${type}
    """
}

process medakaConsensus {
    label "medaka"
    cpus 1
    memory "4 GB"
    input:
        tuple val(meta),
            path("consensus_probs*.hdf"),
            path("align.bam"),
            path("align.bam.bai"),
            path("reference*"),
            val(basecall_model)
    output:
        tuple val(meta), path("${meta.alias}.medaka.fasta.gz")
    shell:
    """
    medaka stitch --threads $task.cpus \
        consensus_probs*.hdf reference* "${meta.alias}.medaka.fasta"

    add_model_to_fasta.sh ${basecall_model} "${meta.alias}.medaka.fasta"
    """
}

process medakaVariant {
    label "medaka"
    cpus 1
    memory "4 GB"
    input:
        tuple val(meta),
            path("consensus_probs*.hdf"),
            path("align.bam"),
            path("align.bam.bai"),
            path("ref.fasta.gz")
    output:
        tuple val(meta), path("${meta.alias}.medaka.vcf.gz"), emit: variants
        tuple val(meta), path("${meta.alias}.variants.stats"), emit: variant_stats
    // note: extension on ref.fasta.gz might not be accurate but shouldn't (?) cause
    //       issues. Also the first step may create an index if not already existing so
    //       the alternative reference.* will break.
    """
    medaka variant ref.fasta.gz consensus_probs*.hdf vanilla.vcf
    bcftools sort vanilla.vcf > vanilla.sorted.vcf

    medaka tools annotate \
        vanilla.sorted.vcf ref.fasta.gz align.bam "${meta.alias}.medaka.vcf"

    bgzip -i "${meta.alias}.medaka.vcf"
    bcftools stats  "${meta.alias}.medaka.vcf.gz" > "${meta.alias}.variants.stats"
    """
}
