process mlstSearch {
    label "mlst"
    cpus 1
    memory "1 GB"
    input:
        tuple val(meta), path("input_genome.fasta.gz")
    output:
        tuple val(meta), path("${meta.alias}.mlst.json")
    script:
    """
    gunzip -c input_genome.fasta.gz > input_genome.fasta
    mlst input_genome.fasta --label ${meta.alias} --json ${meta.alias}.mlst.json
    """
}

process getPointfinderSpecies {
    label "wfbacterialgenomes"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path("${meta.alias}.mlst.json")
    output:
        tuple val(meta), stdout
    script:
    """
    workflow-glue pointfinder_species --mlst_json ${meta.alias}.mlst.json

    """
}

process resfinder {
    label "amr"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta), path("input_genome.fasta.gz"), val(species)
        val resfinder_threshold
        val resfinder_coverage
    output:
        tuple val(meta), path("${meta.alias}_resfinder_results"), val(species)
    script:
    """
    # sed added to remove basecaller config from fasta +  resfinder table 
    gunzip -c input_genome.fasta.gz | sed '/^>/ s/ .*//' > input_genome.fasta

    python -m resfinder \
        -o ${meta.alias}_resfinder_results \
        -j ${meta.alias}_resfinder_results/${meta.alias}_resfinder.json \
        -l ${resfinder_coverage} \
        -t ${resfinder_threshold} \
        --acquired \
        -s "${species}" \
        --point \
        -ifa input_genome.fasta \
        --disinfectant \
        --nanopore
    """
}
process rgi {
	label "rgi"
	cpus 2
	memory "2 GB"
	input :
		tuple val(meta), path("input_genome.fasta.gz")
	output:
		tuple val(meta), path("${meta.alias}_rgi_results.txt")
	script:
	"""
	rgi main -i input_genome.fasta.gz \
	-o ${meta.alias}_rgi_results
	"""
}

process rgi_hamronize {
	label "hamronize"
	cpus 2
	memory "2 GB"
	input :
		tuple val(meta), path("${meta.alias}_rgi_results")
	output:
		tuple val (meta), path("${meta.alias}_hamro_rgi")
	script:
	"""
	hamronize rgi --reference_database_version 6.0.3 --analysis_software_version 6.0.3 --input_file_name ${meta.alias} --output ${meta.alias}_hamro_rgi ${meta.alias}_rgi_results
	"""
}

process amrfinder {
    label "amrfinder"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta), path("input_genome.fasta.gz")
    output:
        tuple val(meta), path("${meta.alias}_amrfinder_results.tsv") 
    script:
    """
    gunzip -c input_genome.fasta.gz > input_genome.fasta
    amrfinder -n input_genome.fasta \
    -o ${meta.alias}_amrfinder_results.tsv
    """
}


process amrfinder_hamronize {
	label "hamronize"
	cpus 2
	memory "2 GB"
	input :
		tuple val(meta), path("${meta.alias}_amrfinder_results.tsv")
	output:
		tuple val (meta), path("${meta.alias}_hamro_amrfinder")
	script:
	"""
	hamronize amrfinderplus --reference_database_version 6.0.3 --analysis_software_version 4.0.3 --input_file_name ${meta.alias} --output ${meta.alias}_hamro_amrfinder ${meta.alias}_amrfinder_results.tsv
	"""
}




process processResfinder {
    // Disinfection not processed yet (CW-2106)
    label "wfbacterialgenomes"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta), path("${meta.alias}_resfinder_results"), val(species)
    output:
        tuple val(meta), path("${meta.alias}.resfinder_results.txt")
    script:
    if (species == "other")
        """
        workflow-glue process_resfinder \
            --resfinder_file ${meta.alias}_resfinder_results/ResFinder_results_tab.txt \
            --output ${meta.alias}.resfinder_results.txt
        """
    else
        """
        workflow-glue process_resfinder \
            --resfinder_file ${meta.alias}_resfinder_results/ResFinder_results_tab.txt \
            --pointfinder_file ${meta.alias}_resfinder_results/PointFinder_results.txt \
            --output ${meta.alias}.resfinder_results.txt \
            --database_location ${meta.alias}_resfinder_results/pointfinder_blast/tmp/
        """
}

process resfinder2 {
    label "amr"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta), path("input_genome.fasta.gz")
        
    output:
        tuple val(meta), path("${meta.alias}_resfinder2.json")
    script:
    """
    # sed added to remove basecaller config from fasta +  resfinder table 
    gunzip -c input_genome.fasta.gz | sed '/^>/ s/ .*//' > ${meta.alias}.fasta

    python -m resfinder \
    	-o ./ -j ${meta.alias}_resfinder2.json \
        --acquired \
        -ifa ${meta.alias}.fasta \
        --disinfectant \
        --nanopore
    """
}

process resfinder_hamronize {
	label "hamronize"
	cpus 2
	memory "2 GB"
	input :
		tuple val(meta), path("${meta.alias}_resfinder2.json")
	output:
		tuple val(meta), path("${meta.alias}_hamro_resfinder")
	script:
	"""
	hamronize resfinder ${meta.alias}_resfinder2.json --output ${meta.alias}_hamro_resfinder
	"""
}

process hamronize {
	label "hamronize"
	cpus 2
	memory "2 GB"
	  
    input:
        tuple val(meta), path("${meta.alias}_hamro_amrfinder")
        tuple val(meta), path("${meta.alias}_hamro_rgi")
        tuple val(meta), path("${meta.alias}_hamro_resfinder")
           
    output:
        tuple val(meta), path("${meta.alias}_hamronized.html")
    
	script:
	"""
	hamronize summarize -o ${meta.alias}_hamronized.html -t interactive ${meta.alias}_hamro_amrfinder ${meta.alias}_hamro_rgi ${meta.alias}_hamro_resfinder	 
	"""
}


process serotyping {
    label "seqsero2"
    cpus 1
    memory "3 GB"
    errorStrategy 'ignore'
    input: 
        tuple val(meta), path("input_genome.fasta.gz"), val(species)
    output:
        tuple val(meta), path("${meta.alias}.serotype_results.tsv")
    script:
    """
    gunzip -c input_genome.fasta.gz > input_genome.fasta

    SeqSero2_package.py \
    -m k \
    -t '4' \
    -i input_genome.fasta \
    -p 1 \
    -b 'mem' \
    -d  output \
    -n ${meta.alias}

    cp -r output/SeqSero_result.tsv "${meta.alias}.serotype_results.tsv"
    """
}   


workflow run_isolates {
   take:
      consensus
      resfinder_threshold
      resfinder_coverage
   main:
        mlst_results = mlstSearch(consensus)
        pointfinder_species = getPointfinderSpecies(mlst_results).map{ meta, species -> [meta, species.trim()] }
        // Added with tuple meta to ensure species tied to correct sample
        resfinder_input = consensus.join(pointfinder_species)
        amr_results = resfinder(resfinder_input, resfinder_threshold, resfinder_coverage)
        processed = processResfinder(amr_results)

        serotype = serotyping(resfinder_input
            | filter { meta, fasta, species -> species == "salmonella" }
        )
        rgi_results = rgi(consensus)
        amrfinder_results = amrfinder(consensus)
        rgi_hamronized = rgi_hamronize(rgi_results)
        amrfinder_hamronized = amrfinder_hamronize(amrfinder_results)
        resfinder_os = resfinder2(consensus)
        resfinder_hamronized = resfinder_hamronize(resfinder_os)
        hamronized = hamronize(amrfinder_hamronized, rgi_hamronized, resfinder_hamronized)
  
   emit:
      amr = amr_results.map{meta, amr, species -> [meta, amr]}
      report_table = processed
      mlst = mlst_results
      serotype = serotype
      rgi = rgi_results
      amrfinder = amrfinder_results
      hamro_rgi = rgi_hamronized
      hamro_amrfinder = amrfinder_hamronized
      hamronize = hamronized 
      
      
      
}
