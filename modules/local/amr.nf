

process resfinderAcquiredOnly {
    label "amr"
    containerOptions '--entrypoint=""'
    input:
        tuple val(meta), path("input_genome.fasta.gz")
    output:
        tuple val(meta), path("${meta.alias}_resfinder_results")
    script:
    """
    gunzip -c input_genome.fasta.gz > input_genome.fasta
    python -m resfinder --acquired -ifa input_genome.fasta --outputPath ${meta.alias}_resfinder_results 
    """
}


process resfinderFull {
    label "amr"
    containerOptions '--entrypoint=""'
    input:
        tuple val(meta), path("input_genome.fasta.gz")
        val species
        val resfinder_threshold
        val resfinder_coverage
    output:
        tuple val(meta), path("${meta.alias}_resfinder_results")
    script: 
        String species_input = species.replace("_", " ");
    """
    gunzip -c -f input_genome.fasta.gz > input_genome.fasta

    python -m resfinder \
        -o ${meta.alias}_resfinder_results \
        -l ${resfinder_coverage} \
        -u \
        -t ${resfinder_threshold} \
        --acquired \
        -s "${species_input}" \
        --point \
        -ifa input_genome.fasta \
        --disinfectant
    """
}


workflow run_amr {
   take:
      consensus
      species
      resfinder_threshold
      resfinder_coverage
   main:
        // If a species does not match the database, the resfinder_full process will fail
        // This can be avoided by using the --ignore_missing_species flag
        // But this is not enabled as otherwise the results may give the wrong impression
        // e.g that point mutations were searched for when they were not
        if (species == "other"){
            amr_results = resfinderAcquiredOnly(consensus)
        } else {
             // if there is a species for the sample then do full amr calling
            amr_results = resfinderFull(consensus, species, resfinder_threshold, resfinder_coverage)
        }
   emit:
      amr = amr_results
}
