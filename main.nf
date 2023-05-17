#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonBuilder

include { fastq_ingress } from './lib/fastqingress'
include { run_amr } from './modules/local/amr'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process readStats {
    label params.process_label
    cpus 1
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
    output:
        path "*readstats.txt", emit: stats
    """
    bamstats align.bam > "${meta.alias}.readstats.txt"
    if [[ \$(wc -l <"${meta.alias}.readstats.txt") -le 1 ]]; then
        echo "No alignments of reads to reference sequence found."
        exit 1
    fi
    """
}


process coverStats {
    label params.process_label
    cpus 2
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
    output:
        path "*fwd.regions.bed.gz", emit: fwd
        path "*rev.regions.bed.gz", emit: rev
        path "*total.regions.bed.gz", emit: all

    """
    mosdepth -n --fast-mode --by 200 --flag 16 -t $task.cpus "${meta.alias}.fwd" align.bam
    mosdepth -n --fast-mode --by 200 --include-flag 16 -t $task.cpus "${meta.alias}.rev" align.bam
    mosdepth -n --fast-mode --by 200 -t $task.cpus "${meta.alias}.total" align.bam
    """
}


process deNovo {
    label params.process_label
    cpus params.threads
    input:
        tuple val(meta), path("reads.fastq.gz")
    output:
        tuple val(meta), path("${meta.alias}.draft_assembly.fasta.gz"), path("${meta.alias}_flye_stats.tsv")
    script:
    """
    flye --nano-raw reads.fastq.gz --out-dir output --threads "${task.cpus}"
    mv output/assembly.fasta "./${meta.alias}.draft_assembly.fasta"
    mv output/assembly_info.txt "./${meta.alias}_flye_stats.tsv"
    bgzip "${meta.alias}.draft_assembly.fasta"
    """
}


process alignReads {
    label params.process_label
    cpus params.threads
    input:
        tuple val(meta), path("reads.fastq.gz"), path("ref.fasta.gz")
    output:
        tuple val(meta), path("*reads2ref.bam"), path("*reads2ref.bam.bai")
    """
    mini_align -i reads.fastq.gz -r ref.fasta.gz -p "${meta.alias}.reads2ref" -t $task.cpus -m
    """
}


process splitRegions {
    // split the bam reference sequences into overlapping sub-regions

    label "medaka"
    cpus 1
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
    output:
        stdout
    """
    #!/usr/bin/env python

    import itertools
    import medaka.common

    regions = itertools.chain.from_iterable(
        x.split(${params.chunk_size}, overlap=1000, fixed_size=False)
        for x in medaka.common.get_bam_regions("align.bam"))
    region_list = []
    for reg in regions:
        # don't ask...just grep &split!
        print("${meta.alias}" + '&split!' + str(reg))
    """
}


// TODO: in a single GPU environment it would be better just
//       to use a single process for the whole bam file. Need
//       to read up on conditional channels

process medakaNetwork {
    // run medaka consensus for each region

    label "medaka"
    cpus 2
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai"), val(reg), val(medaka_model)
    output:
        tuple val(meta), path("*consensus_probs.hdf")
    script:
        def model = medaka_model
    """
    medaka --version
    echo ${model}
    echo ${medaka_model}
    medaka consensus align.bam "${meta.alias}.consensus_probs.hdf" \
        --threads 2 --regions "${reg}" --model ${model}
    """
}


process medakaVariantConsensus {
    // run medaka consensus for each region

    label "medaka"
    cpus 2
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai"), val(reg), val(medaka_model)
    output:
        tuple val(meta), path("*consensus_probs.hdf")
    script:
        def model = medaka_model
    """
    medaka --version
    echo ${model}
    echo ${medaka_model}
    medaka consensus align.bam "${meta.alias}.consensus_probs.hdf" \
        --threads 2 --regions "${reg}" --model ${model}
    """
}


process medakaVariant {
    label "medaka"
    cpus 1
    input:
        tuple val(meta), path("consensus_probs*.hdf"),  path("align.bam"), path("align.bam.bai"), path("ref.fasta.gz")
    output:
        path "${meta.alias}.medaka.vcf.gz", emit: variants
        path "${meta.alias}.variants.stats", emit: variant_stats
    // note: extension on ref.fasta.gz might not be accurate but shouldn't (?) cause issues.
    //       Also the first step may create an index if not already existing so the alternative
    //       reference.* will break
    """
    medaka variant ref.fasta.gz consensus_probs*.hdf vanilla.vcf
    medaka tools annotate vanilla.vcf ref.fasta.gz align.bam "${meta.alias}.medaka.vcf"
    bgzip -i "${meta.alias}.medaka.vcf"
    bcftools stats  "${meta.alias}.medaka.vcf.gz" > "${meta.alias}.variants.stats"
    """
}


process medakaConsensus {
    label "medaka"
    cpus 1
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai"), path("consensus_probs*.hdf"), path("reference*")
    output:
        tuple val(meta), path("${meta.alias}.medaka.fasta.gz")

    """
    medaka stitch --threads $task.cpus consensus_probs*.hdf reference* "${meta.alias}.medaka.fasta"
    bgzip "${meta.alias}.medaka.fasta"
    """
}


process runProkka {
    // run prokka in a basic way on the consensus sequence
    label "prokka"
    cpus params.threads
    input:
        tuple val(meta), path("consensus.fasta.gz")
    output:
        path "*prokka_results/*prokka.gff"

    script:
        def prokka_opts = params.prokka_opts ?: ""
    """
    gunzip -rf consensus.fasta.gz
    prokka $prokka_opts --outdir "${meta.alias}.prokka_results" \
        --cpus $task.cpus --prefix "${meta.alias}.prokka" *consensus.fasta
    """
}


process prokkaVersion {
    label "prokka"
    output:
        path "prokka_version.txt"
    """
    prokka --version | sed 's/ /,/' >> "prokka_version.txt"
    """
}


process medakaVersion {
    label "medaka"
    input:
        path "input_versions.txt"
    output:
        path "medaka_version.txt"
    """
    cat "input_versions.txt" >> "medaka_version.txt"
    medaka --version | sed 's/ /,/' >> "medaka_version.txt"
    """
}


process getVersions {
    label params.process_label
    cpus 1
    input:
        path "input_versions.txt"
    output:
        path "versions.txt"
    """
    cat "input_versions.txt" >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    flye --version | sed 's/^/flye,/' >> versions.txt
    python -c "import pomoxis; print(f'pomoxis,{pomoxis.__version__}')" >> versions.txt
    python -c "import dna_features_viewer; print(f'dna_features_viewer,{dna_features_viewer.__version__}')" >> versions.txt
    """
}


process getParams {
    label params.process_label
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process makeReport {
    label params.process_label
    cpus 1
    input:
        path "versions/*"
        path "params.json"
        path "variants/*"
        val sample_ids
        path "prokka/*"
        path per_read_stats
        path "fwd/*"
        path "rev/*"
        path "total_depth/*"
        path "flye_stats/*"
        path "resfinder/*"
    output:
        path "wf-bacterial-genomes-*.html"
    script:
        report_name = "wf-bacterial-genomes-report.html"
        denovo = params.reference_based_assembly as Boolean ? "" : "--denovo"
        prokka = params.run_prokka as Boolean ? "--prokka" : ""
        resfinder = params.isolates as Boolean ? "--resfinder" : ""
        samples = sample_ids.join(" ")
        String stats_args = \
            (per_read_stats.name == OPTIONAL_FILE.name) ? "" : "--stats $per_read_stats"
    // NOTE: the script assumes the various subdirectories
    """
    workflow-glue report \
    $stats_args \
    $prokka $denovo \
    $resfinder \
    --versions versions \
    --params params.json \
    --output $report_name \
    --sample_ids $samples
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label params.process_label
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}


process lookup_medaka_consensus_model {
    label params.process_label
    input:
        path("lookup_table")
        val basecall_model
    output:
        stdout
    shell:
    '''
    medaka_model=$(workflow-glue resolve_medaka_model lookup_table '!{basecall_model}' "medaka_consensus")
    echo $medaka_model
    '''
}


process lookup_medaka_variant_model {
    label params.process_label
    input:
        path("lookup_table")
        val basecall_model
    output:
        stdout
    shell:
    '''
    medaka_model=$(workflow-glue resolve_medaka_model lookup_table '!{basecall_model}' "medaka_variant")
    echo $medaka_model
    '''
}


// Creates a new directory named after the sample alias and moves the fastcat results
// into it.
process collectFastqIngressResultsInDir {
    label params.process_label
    input:
        // both the fastcat seqs as well as stats might be `OPTIONAL_FILE` --> stage in
        // different sub-directories to avoid name collisions
        tuple val(meta), path(concat_seqs, stageAs: "seqs/*"), path(fastcat_stats,
            stageAs: "stats/*")
    output:
        // use sub-dir to avoid name clashes (in the unlikely event of a sample alias
        // being `seq` or `stats`)
        path "out/*"
    script:
    String outdir = "out/${meta["alias"]}"
    String metaJson = new JsonBuilder(meta).toPrettyString()
    String concat_seqs = \
        (concat_seqs.fileName.name == OPTIONAL_FILE.name) ? "" : concat_seqs
    String fastcat_stats = \
        (fastcat_stats.fileName.name == OPTIONAL_FILE.name) ? "" : fastcat_stats
    """
    mkdir -p $outdir
    echo '$metaJson' > metamap.json
    mv metamap.json $concat_seqs $fastcat_stats $outdir
    """
}

// modular workflow
workflow calling_pipeline {
    take:
        reads
        reference
    main:
        per_read_stats = reads.map {
            it[2] ? it[2].resolve('per-read-stats.tsv') : null
        }
        | collectFile ( keepHeader: true )
        | ifEmpty ( OPTIONAL_FILE )
        input_reads = reads.map { meta, reads, stats -> [meta, reads] }
        sample_ids = reads.map { meta, reads, stats -> meta.alias }
        if (params.reference_based_assembly && !params.reference){
            throw new Exception("Reference based assembly selected, a reference sequence must be provided through the --reference parameter.")
        }
        if (!params.reference_based_assembly){
            log.info("Running Denovo assembly.")
            denovo_assem = deNovo(input_reads)
            named_refs = denovo_assem.map { it -> [it[0], it[1]] }
            read_ref_groups = input_reads.join(named_refs)
        } else {
            log.info("Reference based assembly selected.")
            references = channel.fromPath(params.reference)
            read_ref_groups = input_reads.combine(references)
            named_refs = read_ref_groups.map { it -> [it[0], it[2]] }
        }
        alignments = alignReads(read_ref_groups)
        read_stats = readStats(alignments)
        depth_stats = coverStats(alignments)
        regions = splitRegions(alignments).splitText()
        named_regions = regions.map {
            it -> return tuple(it.split(/&split!/)[0], it.split(/&split!/)[1])
        }

        if(params.medaka_consensus_model) {
            log.warn "Overriding Medaka Consensus model with ${params.medaka_consensus_model}."
            medaka_consensus_model = Channel.fromList([params.medaka_consensus_model])
        }
        else {
            lookup_table = Channel.fromPath("${projectDir}/data/medaka_models.tsv", checkIfExists: true)
            medaka_consensus_model = lookup_medaka_consensus_model(lookup_table, params.basecaller_cfg)
        }
        if(params.medaka_variant_model) {
            log.warn "Overriding Medaka Variant model with ${params.medaka_variant_model}."
            medaka_variant_model = Channel.fromList([params.medaka_variant_model])
        }
        else {
            lookup_table = Channel.fromPath("${projectDir}/data/medaka_models.tsv", checkIfExists: true)
            medaka_variant_model = lookup_medaka_variant_model(lookup_table, params.basecaller_cfg)
        }

        // medaka consensus
        named_alignments = alignments.map{ meta, bam, bai -> [meta.alias, meta, bam, bai] }
        regions_bams = named_alignments.combine(named_regions, by: 0).map{it[1..-1]}
        regions_model = regions_bams.combine(medaka_consensus_model)
        hdfs = medakaNetwork(regions_model)
        hdfs_grouped = alignments.combine(hdfs.groupTuple(), by: 0).join(named_refs)
        consensus = medakaConsensus(hdfs_grouped)
        
        if (!params.reference_based_assembly){
            flye_info = denovo_assem.map { it -> it[2] }
        }else{
            flye_info = Channel.empty()
        }

        // medaka variants
        if (params.reference_based_assembly){
            bam_model = regions_bams.combine(medaka_variant_model)
            hdfs_variant = medakaVariantConsensus(bam_model)
            hdfs_grouped = hdfs_variant.groupTuple().combine(alignments, by: [0]).join(named_refs)
            variant = medakaVariant(hdfs_grouped)
            variants = variant.variant_stats
            vcf_variant = variant.variants
        } else {
            variants = Channel.empty()
            vcf_variant = Channel.empty()
        }

        if (params.run_prokka) {
            prokka = runProkka(consensus)
        } else {
            prokka = Channel.empty()
        }

        // amr calling
        if (params.isolates) {
            run_amr = run_amr(
                consensus,
                params.species,
                "${params.resfinder_threshold}",
                "${params.resfinder_coverage}")
            amr = run_amr.amr 
            amr_results = run_amr.report_table
        } else {
            amr = Channel.empty()
            amr_results = Channel.empty()
        }


        prokka_version = prokkaVersion()
        medaka_version = medakaVersion(prokka_version)
        software_versions = getVersions(medaka_version)
        workflow_params = getParams()

        report = makeReport(
            software_versions.collect(),
            workflow_params,
            variants.collect().ifEmpty(file("${projectDir}/data/OPTIONAL_FILE")),
            sample_ids.collect(),
            prokka.collect().ifEmpty(file("${projectDir}/data/OPTIONAL_FILE")),
            per_read_stats,
            depth_stats.fwd.collect(),
            depth_stats.rev.collect(),
            depth_stats.all.collect(),
            flye_info.collect().ifEmpty(file("${projectDir}/data/OPTIONAL_FILE")),
            amr_results.collect().ifEmpty(file("${projectDir}/data/OPTIONAL_FILE")))
        fastq_stats = reads
        // replace `null` with path to optional file
        | map { [ it[0], it[1] ?: OPTIONAL_FILE, it[2] ?: OPTIONAL_FILE ] }
        | collectFastqIngressResultsInDir
        all_out = variants.concat(
            vcf_variant,
            consensus.map {meta, assembly -> assembly},
            report,
            prokka,
            fastq_stats,
            amr.map {meta, resfinder -> resfinder},
            workflow_params
        )

    emit:
        all_out
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "fastcat_stats": params.wf.fastcat_stats,
        "fastcat_extra_args": ""])

    reference = params.reference
    results = calling_pipeline(samples, reference)
    output(results.all_out)
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
