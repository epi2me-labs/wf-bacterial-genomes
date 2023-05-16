

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


process processResfinderAquired {
    label "wfbacterialgenomes"
    input:
        tuple val(meta), path("${meta.alias}_resfinder_results")
    output:
        path("${meta.alias}.resfinder_results.txt")
    script:
    """
    workflow-glue process_resfinder \
        --resfinder_file ${meta.alias}_resfinder_results/ResFinder_results_tab.txt \
        --output ${meta.alias}.resfinder_results.txt
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


process processResfinderFull {
    label "wfbacterialgenomes"
    input:
        tuple val(meta), path("${meta.alias}_resfinder_results")
    output:
        path("${meta.alias}.resfinder_results.txt")
    script:
    """
    workflow-glue process_resfinder \
        --resfinder_file ${meta.alias}_resfinder_results/ResFinder_results_tab.txt \
        --pointfinder_file ${meta.alias}_resfinder_results/PointFinder_results.txt \
        --output ${meta.alias}.resfinder_results.txt \
        --database_location ${meta.alias}_resfinder_results/pointfinder_blast/tmp/
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
            processed = processResfinderAquired(amr_results)
        } else {
             // if there is a species for the sample then do full amr calling
            amr_results = resfinderFull(consensus, species, resfinder_threshold, resfinder_coverage)
            processed = processResfinderFull(amr_results)
        }
   emit:
      amr = amr_results
      report_table = processed
}
