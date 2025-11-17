process medakaInference {
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
        // TO DO: Add reheader Sam and then update this to use medaka automodel.
        consensus_bact_methyl_compatible_models = [
            'dna_r10.4.1_e8.2_400bps_hac@v4.2.0', 'dna_r10.4.1_e8.2_400bps_sup@v4.2.0', 
            'dna_r10.4.1_e8.2_400bps_hac@v4.3.0', 'dna_r10.4.1_e8.2_400bps_sup@v4.3.0',
            'dna_r10.4.1_e8.2_400bps_hac@v5.0.0', 'dna_r10.4.1_e8.2_400bps_sup@v5.0.0',
            'dna_r10.4.1_e8.2_400bps_hac@v5.2.0', 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0'
        ]
        // Avoid modifying input val, as retry will check the modified val, not original input val
        def model_to_use = basecall_model
        if ((type == "consensus") && (basecall_model in consensus_bact_methyl_compatible_models)){
            model_to_use = "r1041_e82_400bps_bacterial_methylation"       
        }
        else {
            model_to_use = "${basecall_model}:${type}"
        }
    """
    medaka --version
    echo ${basecall_model}
    medaka inference align.bam "${meta.alias}.consensus_probs.hdf" \
        --threads 2 --regions "${region}" --model ${model_to_use}
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
    medaka sequence --threads $task.cpus \
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
    medaka vcf consensus_probs*.hdf ref.fasta.gz vanilla.vcf
    bcftools sort vanilla.vcf > vanilla.sorted.vcf

    medaka tools annotate \
        vanilla.sorted.vcf ref.fasta.gz align.bam "${meta.alias}.medaka.vcf"

    bgzip -i "${meta.alias}.medaka.vcf"
    bcftools stats  "${meta.alias}.medaka.vcf.gz" > "${meta.alias}.variants.stats"
    """
}
