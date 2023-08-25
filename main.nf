#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonBuilder

include { fastq_ingress } from './lib/fastqingress'
include { run_isolates } from './modules/local/isolates'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")
FLYE_MIN_COVERAGE_THRESHOLD = 5

process readStats {
    label "wfbacterialgenomes"
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
    label "wfbacterialgenomes"
    cpus 2
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
    input:
        tuple val(meta), path("reads.fastq.gz")
    output:
        tuple val(meta),
            path("${meta.alias}.draft_assembly.fasta.gz"),
            path("${meta.alias}_flye_stats.tsv"),
            optional: true, emit: asm
        tuple val(meta), env(LOW_COV_FAIL), emit: failed
    script:
    // flye may fail due to low coverage; in this case we don't want to cause the whole
    // workflow to crash --> exit with `0` and don't emit output files
    def flye_opts = params.flye_opts ?: ""
    """
    LOW_COV_FAIL=0
    FLYE_EXIT_CODE=0
    flye $flye_opts --nano-raw reads.fastq.gz --out-dir output --threads "${task.cpus}" || \
    FLYE_EXIT_CODE=\$?

    if [[ \$FLYE_EXIT_CODE -eq 0 ]]; then
        mv output/assembly.fasta "./${meta.alias}.draft_assembly.fasta"
        mv output/assembly_info.txt "./${meta.alias}_flye_stats.tsv"
        bgzip "${meta.alias}.draft_assembly.fasta"
    else
        # flye failed --> check the log to see if low coverage caused the failure
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
            LOW_COV_FAIL=1
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


process medakaVariantHdf {
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
        tuple val(meta), path("${meta.alias}.medaka.vcf.gz"), emit: variants
        tuple val(meta), path("${meta.alias}.variants.stats"), emit: variant_stats
    // note: extension on ref.fasta.gz might not be accurate but shouldn't (?) cause issues.
    //       Also the first step may create an index if not already existing so the alternative
    //       reference.* will break
    """
    medaka variant ref.fasta.gz consensus_probs*.hdf vanilla.vcf
    bcftools sort vanilla.vcf > vanilla.sorted.vcf
    medaka tools annotate vanilla.sorted.vcf ref.fasta.gz align.bam "${meta.alias}.medaka.vcf"
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
    output:
        path "prokka_version.txt"
    """
    prokka --version |& sed 's/ /,/' >> "prokka_version.txt"
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

process mlstVersion {
    label "mlst"
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
    label "wfbacterialgenomes"
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
        path "mlst/*"
    output:
        path "wf-bacterial-genomes-*.html"
    script:
        report_name = "wf-bacterial-genomes-report.html"
        denovo = params.reference_based_assembly as Boolean ? "" : "--denovo"
        prokka = params.run_prokka as Boolean ? "--prokka" : ""
        isolates = params.isolates as Boolean ? "--isolates" : ""
        samples = sample_ids.join(" ")
        String stats_args = \
            (per_read_stats.name == OPTIONAL_FILE.name) ? "" : "--stats $per_read_stats"
    // NOTE: the script assumes the various subdirectories
    """
    workflow-glue report \
    $stats_args \
    $prokka \
    $denovo \
    $isolates \
    --versions versions \
    --params params.json \
    --output $report_name \
    --sample_ids $samples
    """
}


process makePerSampleReports {
    label "wfbacterialgenomes"
    cpus 1
    input:
        path "versions.txt"
        path "params.json"
        tuple val(meta), path("report_files/*")
    output:
        path "${meta.alias}-isolate-report.html"
    script:
        String alias = meta.alias
        String barcode = meta.barcode
        String denovo = params.reference_based_assembly as Boolean ? "" : "--denovo"
        String prokka = params.run_prokka as Boolean ? "--prokka" : ""
        String isolates = params.isolates as Boolean ? "--isolates" : ""
    // the script checks for presence / absence of the various files in `report_files`
    """
    workflow-glue per_sample_report \
        $prokka \
        $denovo \
        $isolates \
        --versions versions.txt \
        --params params.json \
        --output ${meta.alias}-isolate-report.html \
        --sample-alias $alias \
        --sample-barcode $barcode \
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
    label "wfbacterialgenomes"
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
    label "wfbacterialgenomes"
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
    label "wfbacterialgenomes"
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
            deNovo(input_reads)
            // some samples might have failed flye due to low coverage
            deNovo.out.failed.map { meta, failed ->
                if (failed == "1") {
                    log.warn "Flye failed for sample '$meta.alias' due to low coverage."
                }
            }
            named_refs = deNovo.out.asm.map { meta, asm, stats -> [meta, asm] }
            // Nextflow might be run in strict mode (e.g. in CI) which prevents `join`
            // from dropping non-matching entries. We have to use `remainder: true` and
            // filter afterwards instead.
            read_ref_groups = input_reads.join(named_refs, remainder: true).filter {
                meta, reads, asm -> asm
            }
            flye_info = deNovo.out.asm.map { meta, asm, stats -> [meta, stats] }
        } else {
            log.info("Reference based assembly selected.")
            references = channel.fromPath(params.reference)
            read_ref_groups = input_reads.combine(references)
            named_refs = read_ref_groups.map { it -> [it[0], it[2]] }
            flye_info = Channel.empty()
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

        // medaka variants
        if (params.reference_based_assembly){
            bam_model = regions_bams.combine(medaka_variant_model)
            hdfs_variant = medakaVariantHdf(bam_model)
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

        // amr and mlst calling
        if (params.isolates) {
            run_isolates = run_isolates(
                consensus,
                "${params.resfinder_threshold}",
                "${params.resfinder_coverage}")
            mlst = run_isolates.mlst
            amr = run_isolates.amr
            amr_results = run_isolates.report_table
            
            
        } else {
            amr = Channel.empty()
            amr_results = Channel.empty()
            mlst = Channel.empty()
        }

        prokka_version = prokkaVersion()
        medaka_version = medakaVersion(prokka_version)
        mlst_version = mlstVersion(medaka_version)
        software_versions = getVersions(mlst_version)
        workflow_params = getParams()

        report = makeReport(
            software_versions,
            workflow_params,
            variants.map { meta, stats -> stats }.collect().ifEmpty(OPTIONAL_FILE),
            sample_ids.collect(),
            prokka.map{ meta, gff, gbk -> gff }.collect().ifEmpty(OPTIONAL_FILE),
            per_read_stats,
            depth_stats.fwd.map{ meta, depths -> depths }.collect().ifEmpty(OPTIONAL_FILE),
            depth_stats.rev.map{ meta, depths -> depths }.collect().ifEmpty(OPTIONAL_FILE),
            depth_stats.all.map{ meta, depths -> depths }.collect().ifEmpty(OPTIONAL_FILE),
            flye_info.map{ meta, stats -> stats }.collect().ifEmpty(OPTIONAL_FILE),
            amr_results.map{ meta, amr -> amr }.collect().ifEmpty(OPTIONAL_FILE),
            mlst.map{ meta, mlst -> mlst }.collect().ifEmpty(OPTIONAL_FILE))

        if (params.isolates) {
            report_files_per_sample = reads
            | map { meta, reads, stats_dir ->
                [meta, stats_dir.resolve("per-read-stats.tsv")]
            }
            | join(vcf_variant, remainder: true)
            | join(variants, remainder: true)
            | join(prokka, remainder: true)
            | join(depth_stats.fwd, remainder: true)
            | join(depth_stats.rev, remainder: true)
            | join(depth_stats.all, remainder: true)
            | join(flye_info, remainder: true)
            | join(amr, remainder: true)
            | join(mlst, remainder: true)
            | map {
                meta = it[0]
                files = it[1..-1]
                // the empty channels will have resulted in occurrences of `null` in
                // the list produced by the joins --> filter
                [meta, files.findAll { it }]
            }
            perSampleReports = makePerSampleReports(
                software_versions,
                workflow_params,
                report_files_per_sample
            )
        } else {
            perSampleReports = Channel.empty()
        }

        fastq_stats = reads
        // replace `null` with path to optional file
        | map { [ it[0], it[1] ?: OPTIONAL_FILE, it[2] ?: OPTIONAL_FILE ] }
        | collectFastqIngressResultsInDir
        all_out = variants.map{meta, stats -> stats}.concat(
            vcf_variant.map{meta, vcf -> vcf},
            consensus.map {meta, assembly -> assembly},
            report,
            perSampleReports,
            prokka.map{meta, gff, gbk -> [gff, gbk]},
            fastq_stats,
            amr.map {meta, resfinder -> resfinder},
            mlst.map {meta, mlst -> mlst},
            workflow_params,
            software_versions,
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
