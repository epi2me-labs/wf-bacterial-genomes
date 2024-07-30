#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.recursion=true
import groovy.json.JsonBuilder

include { fastq_ingress; xam_ingress  } from './lib/ingress'
include { run_isolates } from './modules/local/isolates'
include {
    medakaHdf as medakaHdf_consensus;
    medakaHdf as medakaHdf_variant;
    medakaConsensus;
    medakaVariant;
} from './modules/local/medaka'

include {
    accumulateCheckpoints;
    ingressCheckpoint;
    assemblyCheckpoint;
    alignmentCheckpoint;
    variantCheckpoint;
    amrCheckpoint;
    annotationCheckpoint;
    perSampleReportingCheckpoint;
    reportingCheckpoint;
} from './modules/local/checkpoints'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")
FLYE_MIN_COVERAGE_THRESHOLD = 5


process readStats {
    label "wfbacterialgenomes"
    cpus 1
    memory "2 GB"
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
    label "wfbacterialgenomes"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
    output:
        tuple val(meta), path("*fwd.regions.bed.gz"), emit: fwd
        tuple val(meta), path("*rev.regions.bed.gz"), emit: rev
        tuple val(meta), path("*total.regions.bed.gz"), emit: all

    """
    mosdepth -n --fast-mode --by 200 --flag 16 -t $task.cpus "${meta.alias}.fwd" align.bam
    mosdepth -n --fast-mode --by 200 --include-flag 16 -t $task.cpus "${meta.alias}.rev" align.bam
    mosdepth -n --fast-mode --by 200 -t $task.cpus "${meta.alias}.total" align.bam
    """
}


process deNovo {
    label "wfbacterialgenomes"
    cpus params.threads
    memory "31 GB"
    input:
        tuple val(meta), path("reads.fastq.gz")
    output:
        tuple val(meta),
            path("${meta.alias}.draft_assembly.fasta.gz"),
            path("${meta.alias}_flye_stats.tsv"),
            optional: true, emit: asm
        tuple val(meta), env(COV_FAIL), emit: failed
    script:
    // flye may fail due to low coverage; in this case we don't want to cause the whole
    // workflow to crash --> exit with `0` and don't emit output files
    def flye_opts = params.flye_opts ?: ""
    def genome_size = params.flye_genome_size ? "--genome-size " + params.flye_genome_size : ""
    def asm_coverage = params.flye_asm_coverage ? "--asm-coverage " + params.flye_asm_coverage : ""
    """
    COV_FAIL=0
    FLYE_EXIT_CODE=0
    flye $flye_opts $genome_size $asm_coverage --nano-hq reads.fastq.gz --out-dir output --threads "${task.cpus}" || \
    FLYE_EXIT_CODE=\$?

    if [[ \$FLYE_EXIT_CODE -eq 0 ]]; then
        mv output/assembly.fasta "./${meta.alias}.draft_assembly.fasta"
        mv output/assembly_info.txt "./${meta.alias}_flye_stats.tsv"
        bgzip "${meta.alias}.draft_assembly.fasta"
    else
        # flye failed --> check the log to check why
        edge_cov=\$(
            grep -oP 'Mean edge coverage: \\K\\d+' output/flye.log \
            || echo $FLYE_MIN_COVERAGE_THRESHOLD
        )
        ovlp_cov=\$(
            grep -oP 'Overlap-based coverage: \\K\\d+' output/flye.log \
            || echo $FLYE_MIN_COVERAGE_THRESHOLD
        )
        if [[
            \$edge_cov -lt $FLYE_MIN_COVERAGE_THRESHOLD ||
            \$ovlp_cov -lt $FLYE_MIN_COVERAGE_THRESHOLD
        ]]; then
            echo -n "Caught Flye failure due to low coverage (either mean edge cov. or "
            echo "overlap-based cov. were below $FLYE_MIN_COVERAGE_THRESHOLD)".
            COV_FAIL=1
        elif grep -q "No disjointigs were assembled" output/flye.log; then
            echo -n "Caught Flye failure due to disjointig assembly."
            COV_FAIL=2
        else
            # exit a subshell with error so that the process fails
            ( exit \$FLYE_EXIT_CODE )
        fi
    fi
    """
}


process alignReads {
    label "wfbacterialgenomes"
    cpus params.threads
    memory "8 GB"
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
    memory "4 GB"
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
    output:
        path "output.txt"
    """
    #!/usr/bin/env python

    import itertools
    import medaka.common

    regions = itertools.chain.from_iterable(
        x.split(${params.chunk_size}, overlap=1000, fixed_size=False)
        for x in medaka.common.get_bam_regions("align.bam"))
    region_list = []
    with open("output.txt", "w") as outfile:
        for reg in regions:
            # don't ask...just grep &split!
            outfile.write("${meta.alias}" + '&split!' + str(reg) + "\\n")
    """
}

process runProkka {
    // run prokka in a basic way on the consensus sequence
    label "prokka"
    cpus params.threads
    memory "4 GB"
    input:
        tuple val(meta), path("consensus.fasta.gz")
    output:
        tuple val(meta), path("*prokka_results/*prokka.gff"), path("*prokka_results/*prokka.gbk")

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
    cpus 1
    memory "2 GB"
    output:
        path "prokka_version.txt"
    """
    prokka --version |& sed 's/ /,/' >> "prokka_version.txt"
    """
}


process medakaVersion {
    label "medaka"
    cpus 1
    memory "2 GB"
    input:
        path "input_versions.txt"
    output:
        path "medaka_version.txt"
    """
    cat "input_versions.txt" >> "medaka_version.txt"
    medaka --version | sed 's/ /,/' >> "medaka_version.txt"
    """
}

process mlstVersion {
    label "mlst"
    cpus 1
    memory "2 GB"
    input:
        path "input_version.txt"
    output:
        path "mlst_version.txt"
    """
    cat "input_version.txt" >> "mlst_version.txt"
    mlst --version | sed 's/ /,/' >> "mlst_version.txt"
    """
}



process getVersions {
    label "wfbacterialgenomes"
    cpus 1
    memory "2 GB"
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
    """
}


process getParams {
    label "wfbacterialgenomes"
    cpus 1
    memory "2 GB"
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process collect_results {
    label "wfbacterialgenomes"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path("report_files/*")
        path("params.json")
    output:
        path "${meta.alias}.json"
    script:
        String alias = meta.alias
        String barcode = meta.barcode
        String type = meta.type
    """
    workflow-glue collect_results \
        --output ${alias}.json \
        --alias $alias \
        --barcode $barcode \
        --params params.json \
        --type $type \
        --data_dir report_files
    """
}


process createRunModel {
    label "wfbacterialgenomes"
    cpus 1
    memory "15 GB"
    input:
        path "sample_results/*"
        val metadata
    output:
        path "results.json"
    script:
    metaJson = new JsonBuilder(metadata).toString()
    """
    workflow-glue create_run_model \
        --jsons sample_results/* \
        --metadata '${metaJson}' \
        --output results.json
    """
}


process makeReport {
    label "wfbacterialgenomes"
    cpus 1
    memory "15 GB"
    input:
        path "versions/*"
        path "params.json"
        path "variants/*"
        val sample_ids
        path "prokka/*"
        val sample_ids_with_stats
        path "per_read_stats/?.gz"
        path "fwd/*"
        path "rev/*"
        path "total_depth/*"
        path "flye_stats/*"
        path "resfinder/*"
        path "mlst/*"
        path "serotype/*"
        path client_fields
    output:
        path "wf-bacterial-genomes-*.html"
    script:
        report_name = "wf-bacterial-genomes-report.html"
        denovo = params.reference_based_assembly as Boolean ? "" : "--denovo"
        prokka = params.run_prokka as Boolean ? "--prokka" : ""
        isolates = params.isolates as Boolean ? "--isolates" : ""
        samples = sample_ids.join(" ")
        sample_ids_with_stats_arg = sample_ids_with_stats ? \
            "--sample_ids_with_stats ${sample_ids_with_stats.join(" ")}" : ""
        client_fields_args = client_fields.name == OPTIONAL_FILE.name ? "" : "--client_fields $client_fields"
    // NOTE: the script assumes the various subdirectories
    """
    workflow-glue report \
    --stats per_read_stats/* \
    $sample_ids_with_stats_arg \
    $prokka \
    $denovo \
    $isolates \
    --versions versions \
    --params params.json \
    --output $report_name \
    --sample_ids $samples \
   $client_fields_args 
    """
}


process makePerSampleReports {
    label "wfbacterialgenomes"
    cpus 1
    memory "15 GB"
    input:
        path "versions.txt"
        path "params.json"
        tuple val(meta), path("report_files/*")
    output:
        tuple val(meta), path("${meta.alias}-isolate-report.html")
    script:
        String barcode = meta.barcode
        String denovo = params.reference_based_assembly as Boolean ? "" : "--denovo"
    // the script checks for presence / absence of the various files in `report_files`
    """
    workflow-glue per_sample_report \
        $denovo \
        --versions versions.txt \
        --params params.json \
        --output ${meta.alias}-isolate-report.html \
        --sample-alias ${meta.alias} \
        --sample-barcode $barcode \
        --data_dir report_files \
        --wf-session $workflow.sessionId \
        --wf-version $workflow.manifest.version
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wfbacterialgenomes"
    cpus 1
    memory "2 GB"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}

// Creates a new directory named after the sample alias and moves the fastcat results
// into it.
process collectFastqIngressResultsInDir {
    label "wfbacterialgenomes"
    cpus 1
    memory "2 GB"
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
        reads.branch { meta, reads, stats -> 
            reads : reads != null
                return [ meta, reads ]
            no_reads : reads == null
                return [ meta, OPTIONAL_FILE ]
        }.set{input_reads}
        
        ingress_checkpoint = ingressCheckpoint(
            input_reads.reads | map { meta, reads -> [ meta, "complete" ] }
            | mix (input_reads.no_reads | map { meta, reads -> [ meta, "not-met" ] } )
        )
        
        sample_ids = reads.map { meta, reads, stats -> meta.alias }
        metadata = reads.map { meta, reads, stats -> meta } | toList()
        definitions = projectDir.resolve("./output_definition.json").toString()
        client_fields = params.client_fields && file(params.client_fields).exists() ? file(params.client_fields) : OPTIONAL_FILE

        // get basecall models: we use `params.override_basecaller_cfg` if present;
        // otherwise use `meta.basecall_models[0]` (there should only be one value in
        // the list because we're running ingress with `allow_multiple_basecall_models:
        // false`; note that `[0]` on an empty list returns `null`)
        basecall_models = reads.map { meta, reads, stats ->
            String basecall_model = \
                params.override_basecaller_cfg ?: meta.basecall_models[0]
            if (!basecall_model) {
                error "Found no basecall model information in the input data for " + \
                    "sample '$meta.alias'. Please provide it with the " + \
                    "`--override_basecaller_cfg` parameter."
            }
            [meta, basecall_model]
        }

        if (params.reference_based_assembly && !params.reference){
            throw new Exception("Reference based assembly selected, a reference sequence must be provided through the --reference parameter.")
        }
        if (!params.reference_based_assembly){
            log.info("Running Denovo assembly.")
            deNovo(input_reads.reads)
            // some samples might have failed flye due to low coverage
            deNovo.out.failed.map { meta, failed ->
                if (failed == "1") {
                    log.warn "Flye failed for sample '$meta.alias' due to low coverage."
                } else if (failed == "2"){
                    log.warn "Flye failed for sample '$meta.alias' as no disjointigs were assembled."
                }
            }

            // Creat channel of failed samples for checkpoints "not-met"
            failed_samples = input_reads.no_reads.mix(
                deNovo.out.failed | filter { meta, failed -> failed != "0"}
            ) | map { meta, field -> [ meta, "not-met" ] }
            named_refs = deNovo.out.asm.map { meta, asm, stats -> [meta, asm] }
            // Nextflow might be run in strict mode (e.g. in CI) which prevents `join`
            // from dropping non-matching entries. We have to use `remainder: true` and
            // filter afterwards instead.
            read_ref_groups = input_reads.reads.join(named_refs, remainder: true).filter {
                meta, reads, asm -> asm
            }
            flye_info = deNovo.out.asm.map { meta, asm, stats -> [meta, stats] }
        } else {
            log.info("Reference based assembly selected.")
            references = channel.fromPath(params.reference)
            read_ref_groups = input_reads.reads.combine(references)
            named_refs = read_ref_groups.map { it -> [it[0], it[2]] }
            flye_info = Channel.empty()
            failed_samples = input_reads.no_reads 
                | map { meta, reads -> [ meta, "not-met" ] }
        }

        alignments = alignReads(read_ref_groups)
        
        // Checkpoint 1 - Alignment
        alignment_checkpoint = alignmentCheckpoint(alignments
        | concat( failed_samples
        | map {meta, status -> [ meta, OPTIONAL_FILE, OPTIONAL_FILE ] } ) )

        read_stats = readStats(alignments)
        depth_stats = coverStats(alignments)
        regions = splitRegions(alignments).splitText()
        named_regions = regions.map {
            it -> return tuple(it.split(/&split!/)[0], it.split(/&split!/)[1])
        }

        // medaka consensus
        named_alignments = alignments.map{ meta, bam, bai -> [meta.alias, meta, bam, bai] }
        // use `sample_id` to combine here
        regions_bams = named_alignments.combine(named_regions, by: 0).map{it[1..-1]}
        regions_model = regions_bams.combine(basecall_models, by: 0)
        // the `.combine`s below use the meta map (and not sample id)
        consensus = medakaHdf_consensus(regions_model, "consensus")
        | groupTuple
        | combine(alignments, by: 0)
        | combine(named_refs, by: 0)
        | combine(basecall_models, by: 0)
        | medakaConsensus

        // Checkpoint 2 - Assembly
        assembly_checkpoint = assemblyCheckpoint(consensus
        | concat (failed_samples
        | map { meta, status -> [ meta, OPTIONAL_FILE ] } ))

        // medaka variants
        if (params.reference_based_assembly){
            medakaHdf_variant(regions_model, "variant")
            | groupTuple
            | combine(alignments, by: 0)
            | combine(named_refs, by: 0)
            | medakaVariant

            vcf_stats = medakaVariant.out.variant_stats
            vcf_variant = medakaVariant.out.variants
            vcf_status = vcf_variant
                | map { meta, variants -> [ meta, "complete" ] }

        } else {
            vcf_stats = Channel.empty()
            vcf_variant = Channel.empty() 
            vcf_status = reads
                | map { meta, reads , stats -> [ meta, "not-met" ] }
        }

        // Checkpoint 3 - variants
        variant_checkpoint = variantCheckpoint(vcf_status
        | mix( failed_samples )
        | unique() )


        if (params.run_prokka) {
            prokka = runProkka(consensus)
            prokka_status = prokka |
                map { meta, gff, gbk  -> [ meta, "complete" ] }
        } else {
            prokka = Channel.empty()
            prokka_status = reads |
                map { meta, reads, stats -> [ meta, "not-met" ] }
        }

        // Checkpoint 4 - annotations
        annotation_checkpoint = annotationCheckpoint(prokka_status
        | mix( failed_samples )
        | unique() )

        // amr and mlst calling
        if (params.isolates) {
            run_isolates = run_isolates(
                consensus,
                "${params.resfinder_threshold}",
                "${params.resfinder_coverage}")
            mlst = run_isolates.mlst
            amr = run_isolates.amr
            amr_results = run_isolates.report_table
            serotype = run_isolates.serotype
            amr_status = amr_results |
                map { meta, resfinder -> [ meta, "complete" ] }
            
        } else {
            amr = Channel.empty()
            amr_results = Channel.empty()
            mlst = Channel.empty()
            serotype = Channel.empty()
            amr_status = reads |
                map { meta, reads, stats -> [ meta, "not-met" ] }
        }

        // Checkpoint 5 - AMR / isolates
        amr_checkpoint = amrCheckpoint(amr_status
        | mix( failed_samples )
        | unique() )

        prokka_version = prokkaVersion()
        medaka_version = medakaVersion(prokka_version)
        mlst_version = mlstVersion(medaka_version)
        software_versions = getVersions(mlst_version)
        workflow_params = getParams()

        // Taken from per sample reports to fill in wf.Sample
        // This is a temporary solution before reporting is done with results.json CW-3217
        report_files_per_sample = reads | filter {meta, reads, stats -> reads != null }
            | map { meta, reads, stats_dir -> 
                
                [meta, stats_dir ? stats_dir : null]
            }
            | join(vcf_variant, remainder: true)
            | join(vcf_stats, remainder: true)
            | join(prokka, remainder: true)
            | join(depth_stats.fwd, remainder: true)
            | join(depth_stats.rev, remainder: true)
            | join(depth_stats.all, remainder: true)
            | join(flye_info, remainder: true)
            | join(amr, remainder: true)
            | join(mlst, remainder: true)
            | join(serotype, remainder: true)
            | map {
                meta = it[0]
                files = it[1..-1]
                // the empty channels will have resulted in occurrences of `null` in
                // the list produced by the joins --> filter
                [meta, files.findAll { it }]
            }
        
        sample_jsons = collect_results(report_files_per_sample, workflow_params)

        run_model = createRunModel(
            sample_jsons.collect(),
            metadata
        )

        // get the samples that have stats (to keep the aliases and stats in the same
        // order for the report)
        samples_with_stats = reads.map { meta, reads, stats_dir ->
            if (stats_dir) {
                [meta.alias, stats_dir.resolve('per-read-stats.tsv.gz')]
            }
        }
        | multiMap { alias, stats ->
            alias: alias
            stats: stats
        }

        report = makeReport(
            software_versions,
            workflow_params,
            vcf_stats.map { meta, stats -> stats }.collect().ifEmpty(OPTIONAL_FILE),
            sample_ids.collect(),
            prokka.map{ meta, gff, gbk -> gff }.collect().ifEmpty(OPTIONAL_FILE),
            samples_with_stats.alias.toList(),
            samples_with_stats.stats.toList(),
            depth_stats.fwd.map{ meta, depths -> depths }.collect().ifEmpty(OPTIONAL_FILE),
            depth_stats.rev.map{ meta, depths -> depths }.collect().ifEmpty(OPTIONAL_FILE),
            depth_stats.all.map{ meta, depths -> depths }.collect().ifEmpty(OPTIONAL_FILE),
            flye_info.map{ meta, stats -> stats }.collect().ifEmpty(OPTIONAL_FILE),
            amr_results.map{ meta, amr -> amr }.collect().ifEmpty(OPTIONAL_FILE),
            mlst.map{ meta, mlst -> mlst }.collect().ifEmpty(OPTIONAL_FILE),
            serotype.map{ meta, sero -> sero}.collect().ifEmpty(OPTIONAL_FILE),
            client_fields)
        
        // Checkpoint 6 - report
        reporting_checkpoint = reportingCheckpoint(report)


        if (params.isolates) {
            perSampleReports = makePerSampleReports(
                software_versions,
                workflow_params,
                report_files_per_sample
            )
            per_sample_report_status = perSampleReports 
                | map { meta, report -> [ meta, "complete" ] }
        } else {
            perSampleReports = Channel.empty()
            per_sample_report_status = reads |
                map { meta, reads, stats -> [ meta, "not-met" ] }
        }
        
        // Checkpoint 7 - per sample report
        per_sample_reporting_checkpoint = perSampleReportingCheckpoint(per_sample_report_status
        | mix( failed_samples )
        | unique() )

        accumulateCheckpoints.scan(
        ingress_checkpoint.mix(
            alignment_checkpoint,
            assembly_checkpoint,
            variant_checkpoint,
            annotation_checkpoint,
            amr_checkpoint,
            reporting_checkpoint,
            per_sample_reporting_checkpoint
        ),
        metadata,
        definitions
    )

        fastq_stats = reads
        // replace `null` with path to optional file
        | map { [ it[0], it[1] ?: OPTIONAL_FILE, it[2] ?: OPTIONAL_FILE ] }
        | collectFastqIngressResultsInDir
        all_out = vcf_stats.map{meta, stats -> stats}.concat(
            vcf_variant.map {meta, vcf -> vcf},
            consensus.map {meta, assembly -> assembly},
            report,
            perSampleReports.map {meta, report -> report},
            prokka.map{meta, gff, gbk -> [gff, gbk]},
            fastq_stats,
            amr.map {meta, resfinder -> resfinder},
            mlst.map {meta, mlst -> mlst},
            workflow_params,
            software_versions,
            run_model,
            serotype.map { meta, sero -> sero }
        )

    emit:
        all_out
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

      File checkpoints_file = new File("checkpoints.json");  

    if (checkpoints_file.exists() == true && workflow.resume == false){
      checkpoints_file.delete()
    }

    // warn the user if overriding the basecall models found in the inputs
    if (params.override_basecaller_cfg) {
        log.warn \
            "Overriding basecall model with '${params.override_basecaller_cfg}'."
    }

    String fastcat_extra_args = params.min_read_length ? " -a $params.min_read_length " : ""

    Map ingress_args = [
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "stats":true,
        "per_read_stats":true,
        "fastcat_extra_args":fastcat_extra_args,
        "allow_multiple_basecall_models": false,
    ]

    if (params.fastq){
        samples = fastq_ingress(ingress_args + [
            "input":params.fastq,
        ])
    } else {
        samples = xam_ingress(ingress_args + [
            "input":params.bam,
            "return_fastq":true,
        ])
    } 

   
    reference = params.reference
    results = calling_pipeline(samples, reference)

    results.all_out
    | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}

