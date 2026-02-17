#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.recursion=true
import groovy.json.JsonBuilder

include { fastq_ingress; xam_ingress  } from './lib/ingress'
include { run_isolates } from './modules/local/isolates'
include {
    medakaInference as medakaInference_consensus;
    medakaInference as medakaInference_variant;
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
    plasmidIdCheckpoint;
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
    memory { task.attempt == 1 ? "31 GB" : task.attempt == 2 ? "63 GB" : "127 GB" }
    maxRetries 3
    errorStrategy 'retry'
    input:
        tuple val(meta), path("reads.fastq.gz")
    output:
        tuple val(meta),
            path("${meta.alias}.draft_assembly.fasta.gz"),
            path("${meta.alias}.flye_stats.tsv"),
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
        mv output/assembly_info.txt "./${meta.alias}.flye_stats.tsv"
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
        tuple val(meta), path("*.bam"), path("*.bam.bai")
    """
    mini_align -i reads.fastq.gz -r ref.fasta.gz -p "${meta.alias}" -t $task.cpus -m
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

process runDnaapler {
    label "wfbacterialgenomes"
    cpus params.threads
    memory "4 GB"
    input:
        tuple val(meta), path("assembly.fasta.gz")
    output:
        tuple val(meta), path("dnaapler_output/${meta.alias}.draft_assembly.fasta.gz"), emit: assembly
        tuple val(meta), env(DNAAPLER_EXIT_CODE), env(DNAAPLER_STDERR), emit: exit_status
    script:
    """
    gunzip -c assembly.fasta.gz > assembly.fasta
    set +e
    dnaapler all \
        --input assembly.fasta \
        --output dnaapler_output \
        --prefix ${meta.alias}_dnaapler \
        --threads $task.cpus 2> stderr.log
    DNAAPLER_EXIT_CODE=\$?
    set -e
    if [[ \$DNAAPLER_EXIT_CODE -ne 0 ]]; then
        DNAAPLER_STDERR=\$(cat stderr.log 2>/dev/null || echo "No stderr available")
    fi
    if [[ \$DNAAPLER_EXIT_CODE -eq 0 ]] && [[ -s "dnaapler_output/${meta.alias}_dnaapler_reoriented.fasta" ]]; then
        # If exit code 0 and non-empty output
        bgzip -c "dnaapler_output/${meta.alias}_dnaapler_reoriented.fasta" > "dnaapler_output/${meta.alias}.draft_assembly.fasta.gz"
    else
        # On dnaapler failure, output original .fasta
        bgzip -c assembly.fasta > "dnaapler_output/${meta.alias}.draft_assembly.fasta.gz"
    fi
    """
}


process runBakta {
    // run Bakta on the consensus sequence
    label "bakta"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta), path("consensus.fasta.gz")
        path bakta_db
    output:
        tuple val(meta), 
            path("${meta.alias}.bakta_results/${meta.alias}.bakta.gff3"), 
            path("${meta.alias}.bakta_results/${meta.alias}.bakta.gbff"),
            optional: true, emit: annot
        tuple val(meta), 
            env(BAKTA_EXIT_CODE), 
            env(STDERR),
            emit: exit_status
    script:
    def bakta_opts = params.bakta_opts ?: ""
    def db_arg = bakta_db.name != 'OPTIONAL_FILE' ? "--db ${bakta_db}" : ""
    """
    gunzip -c consensus.fasta.gz > consensus.fasta
    set +e
    bakta ${db_arg} \
        ${bakta_opts} \
        --keep-contig-headers \
        --output "${meta.alias}.bakta_results" \
        --threads $task.cpus \
        --prefix "${meta.alias}.bakta" \
        --skip-plot \
        *consensus.fasta 2> stderr.log
    BAKTA_EXIT_CODE=\$?
    set -e
    if [[ \$BAKTA_EXIT_CODE -ne 0 ]]; then
        STDERR=\$(cat stderr.log 2>/dev/null || echo "No stderr available")
    else
        STDERR=""
    fi
    """
}


process runMobSuite {
    label "mobsuite"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta), path("assembly.fasta.gz")
    output:
        tuple val(meta), 
            path("${meta.alias}_mob_results/"),
            optional: true, emit: results
        tuple val(meta), 
            env(MOB_EXIT_CODE), 
            env(MOB_STDERR),
            emit: exit_status
    script:
    def mob_opts = params.plasmid_id_opts ?: ""
    """
    gunzip -c assembly.fasta.gz > assembly.fasta
    mkdir -p "${meta.alias}_mob_results"
    
    set +e
    mob_recon \
        ${mob_opts} \
        --infile assembly.fasta \
        --outdir "${meta.alias}_mob_results" \
        --force \
        --num_threads ${task.cpus} 2> stderr.log
    MOB_EXIT_CODE=\$?
    set -e
    
    if [[ \$MOB_EXIT_CODE -ne 0 ]]; then
        MOB_STDERR=\$(cat stderr.log 2>/dev/null || echo "No stderr available")
    else
        MOB_STDERR=""
    fi
    """
}

process baktaVersion {
    label "bakta"
    cpus 1
    memory "2 GB"
    output:
        path "bakta_version.txt"
    """
    bakta --version | sed 's/ /,/' >> "bakta_version.txt"
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
    bcftools --version | head -n 1 | sed 's/ /,/' >> "medaka_version.txt"
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

process sourmashVersion {
    label "sourmash"
    cpus 1
    memory "2 GB"
    input:
        path "input_version.txt"
    output:
        path "sourmash_version.txt"
    """
    cat "input_version.txt" >> "sourmash_version.txt"
    sourmash --version | sed 's/ /,/' >> "sourmash_version.txt"
    """
}

process mobsuiteVersion {
    label "mobsuite"
    cpus 1
    memory "2 GB"
    input:
        path "input_version.txt"
    output:
        path "mobsuite_version.txt"
    """
    cat "input_version.txt" >> "mobsuite_version.txt"
    mob_recon --version | sed 's/ /,/' >> "mobsuite_version.txt"
    """
}

process resfinderVersion {
    label "amr"
    cpus 1
    memory "2 GB"
    input:
        path "input_version.txt"
    output:
        path "resfinder_version.txt"
    """
    cat "input_version.txt" >> "resfinder_version.txt"
    python -m resfinder --version 2>&1 | sed 's/^/resfinder,/' >> "resfinder_version.txt"
    """
}

process seqseroVersion {
    label "seqsero2"
    cpus 1
    memory "2 GB"
    input:
        path "input_version.txt"
    output:
        path "seqsero_version.txt"
    """
    cat "input_version.txt" >> "seqsero_version.txt"
    SeqSero2_package.py --version 2>&1 | sed 's/SeqSero2_package.py/SeqSero2/' | sed 's/ /,/' >> "seqsero_version.txt"
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
    bamstats --version | sed 's/^/bamstats,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    flye --version | sed 's/^/flye,/' >> versions.txt
    dnaapler --version | head -1 | sed 's/, version /,/' >> versions.txt
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
        tuple val(meta), path("${meta.alias}.json")
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
    label "wf_common"
    cpus 1
    memory "15 GB"
    input:
        path "versions/*"
        path "params.json"
        path "variants/*"
        val sample_ids
        val bakta_enabled 
        path "per_read_stats/*"
        path "fwd/*"
        path "rev/*"
        path "total_depth/*"
        path results
        path client_fields
    output:
        path "wf-bacterial-genomes-*.html"
    script:
        report_name = "wf-bacterial-genomes-report.html"
        denovo = params.reference_based_assembly as Boolean ? "" : "--denovo"
        bakta_arg = bakta_enabled && params.run_bakta ? "--bakta" : ""
        isolates = params.isolates as Boolean ? "--isolates" : ""
        plasmid_id = params.run_plasmid_id as Boolean ? "--plasmid_id" : ""  
        samples = sample_ids.join(" ")
        client_fields_args = client_fields.name == OPTIONAL_FILE.name ? "" : "--client_fields $client_fields"
    // NOTE: the script assumes the various subdirectories
    """
    workflow-glue report \
    $bakta_arg \
    $denovo \
    $isolates \
    $plasmid_id \
    --versions versions \
    --params params.json \
    --output $report_name \
    --sample_ids $samples \
    --results $results \
    $client_fields_args \
    --wf_version ${workflow.manifest.version}
    """
}


process makePerSampleReports {
    label "wf_common"
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
        --sample_alias ${meta.alias} \
        --sample_barcode $barcode \
        --data_dir report_files \
        --wf_session $workflow.sessionId \
        --wf_version $workflow.manifest.version
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
            reads : meta.n_seqs > 0
                return [ meta, reads ]
            no_reads : meta.n_seqs == null || meta.n_seqs == 0
                return [ meta, OPTIONAL_FILE ]
        }.set{input_reads}
        
        ingress_checkpoint = ingressCheckpoint(
            input_reads.reads | map { meta, reads -> [ meta, "complete" ] }
            | mix (input_reads.no_reads | map { meta, reads -> [ meta, "not-met" ] } )
        )


        // get basecall models: we use `params.override_basecaller_cfg` if present;
        // otherwise use `meta.basecall_models[0]` (there should only be one value in
        // the list because we're running ingress with `allow_multiple_basecall_models:
        // false`; note that `[0]` on an empty list returns `null`)

        basecall_models_initial = input_reads.reads.map { meta, reads ->
            String basecall_model = \
                params.override_basecaller_cfg ?: meta.basecall_models[0]
            if (!basecall_model) {
                error "Found no basecall model information in the input data for " + \
                    "sample '$meta.alias'. Please provide it with the " + \
                    "`--override_basecaller_cfg` parameter."
            }
            [meta, basecall_model]
        }

        sample_ids = reads.map { meta, reads, stats -> meta.alias }
        metadata = reads.map { meta, reads, stats -> meta } | toList()
        definitions = projectDir.resolve("./output_definition.json").toString()
        client_fields = params.client_fields && file(params.client_fields).exists() ? file(params.client_fields) : OPTIONAL_FILE

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
            
            // Attempt to reorient successful Flye assemblies to dnaA/repA
            flye_assemblies = deNovo.out.asm.map { meta, asm, stats -> [meta, asm] }
            dnaapler_output = runDnaapler(flye_assemblies)
            dnaapler_output.exit_status.subscribe { meta, exit_code, stderr ->
                if (exit_code != "0") {
                    log.warn "Dnaapler failed for sample '${meta.alias}': ${stderr}"
                }
            }
            named_refs = dnaapler_output.assembly

            // Nextflow might be run in strict mode (e.g. in CI) which prevents `join`
            // from dropping non-matching entries. We have to use `remainder: true` and
            // filter afterwards instead.
            read_ref_groups = input_reads.reads.join(named_refs, remainder: true).filter {
                meta, reads, asm -> asm
            }
            flye_info = deNovo.out.asm.map { meta, asm, stats -> [meta, stats] }
        } else {
            log.info("Reference based assembly selected.")
            references = Channel.fromPath(params.reference)
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
        

        // Filter out samples that failed assembly
        // Join will create [meta, basecall, null] for passed samples
        // and [meta, null, "not-met"] for failed samples
        basecall_models = basecall_models_initial
        | join( failed_samples, failOnMismatch: false, remainder: true)
        | filter{ meta, basecall, failed -> failed == null }
        | map { meta, basecall, failed -> [meta, basecall] }

        // medaka consensus
        named_alignments = alignments.map{ meta, bam, bai -> [meta.alias, meta, bam, bai] }
        // use `sample_id` to combine here
        regions_bams = named_alignments.combine(named_regions, by: 0).map{it[1..-1]}
        regions_model = regions_bams.combine(basecall_models, by: 0)
        // the `.combine`s below use the meta map (and not sample id)
        consensus= medakaInference_consensus(regions_model, "consensus")
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
            medakaInference_variant(regions_model, "variant")
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

        if (params.run_plasmid_id) {
            mob_output = runMobSuite(consensus)
            
            // Log failed runs
            mob_output.exit_status.subscribe { meta, exit_code, stderr ->
                if (exit_code != "0") {
                    log.warn "Mob-suite failed for sample '${meta.alias}': ${stderr}"
                }
            }
            
            // Filter successful runs
            mobsuite_results = mob_output.exit_status
                .join(mob_output.results, by: 0, remainder: true, failOnMismatch: false)
                .filter { meta, exit_code, stderr, results = null ->
                    def valid_output = results?.exists()
                    return exit_code == "0" && valid_output
                }
                .map { meta, exit_code, stderr, results -> [meta, results] }
            
            mobsuite_status = mobsuite_results
                .map { meta, results -> [meta, "complete"] }
                .mix(failed_samples.map { meta, status -> [meta, "not-met"] })
                .unique()
        } else {
            mobsuite_results = Channel.empty()
            mobsuite_status = reads | map { meta, reads, stats -> [meta, "not-met"] }
        }

        // Checkpoint 4 - plasmid identification
        plasmid_checkpoint = plasmidIdCheckpoint(mobsuite_status
            | mix(failed_samples)
            | unique())

        if (params.run_bakta) {
            bakta_db = params.bakta_db ? 
                Channel.fromPath(params.bakta_db, checkIfExists: true) : 
                Channel.value(OPTIONAL_FILE)

            bakta_output = runBakta(consensus, bakta_db)
            // Check if bakta annotation succedds for each sample and log failed samples with error message
            bakta_output.exit_status.subscribe { meta, exit_code, stderr_tail ->
                if (exit_code != "0") {
                    log.warn "Bakta annotation failed for sample '${meta.alias}': ${stderr_tail}"
                }
            }        
            // Pass on successful runs only -> exit code 0 + files exist
            bakta = bakta_output.exit_status
                .join(bakta_output.annot, by: 0, remainder: true, failOnMismatch: false)
                .filter { meta, exit_code, stderr_tail, gff3 = null, gbff = null ->
                    // Check that output files exist and are not empty
                    def valid_output = gff3?.size() > 0 && gbff?.size() > 0
                    return exit_code == "0" && valid_output
                }
                .map { meta, exit_code, stderr_tail, gff3, gbff -> [meta, gff3, gbff] }
            
            bakta_status = bakta
                .map { meta, gff3, gbff -> [meta, "complete"] }
                .mix(failed_samples.map { meta, status -> [meta, "not-met"] })
                .unique()
        }
        else {
            bakta = Channel.empty()
            bakta_status = reads |
                map { meta, reads, stats -> [ meta, "not-met" ] }
        }
        bakta_files = bakta.map{meta, gff3, gbff -> gff3}.collect().ifEmpty(OPTIONAL_FILE)

        // Checkpoint 5 - annotations
        annotation_checkpoint = annotationCheckpoint(bakta_status
        | mix( failed_samples )
        | unique() )

        resfinder_db = params.resfinder_db ? 
            Channel.fromPath(params.resfinder_db, checkIfExists: true) : 
            Channel.value(OPTIONAL_FILE)
        pointfinder_db = params.pointfinder_db ? 
            Channel.fromPath(params.pointfinder_db, checkIfExists: true) : 
            Channel.value(OPTIONAL_FILE)
        db_mapping_file = params._sourmash_pointfinder_mapping ? 
            Channel.fromPath(params._sourmash_pointfinder_mapping, checkIfExists: true) : 
            Channel.value(OPTIONAL_FILE)

        // amr and mlst calling
        if (params.isolates) {
            run_isolates = run_isolates(
                consensus,
                params._sourmash_db_exclude_list,
                "${params.resfinder_threshold}",
                "${params.resfinder_coverage}",
                params.pointfinder_ignore_indels,
                params.pointfinder_ignore_stop_codons,
                resfinder_db,
                pointfinder_db,
                db_mapping_file)
            mlst = run_isolates.mlst
            taxonomy_results = run_isolates.taxonomy
            sourmash_exclude = run_isolates.sourmash_excluded_genomes
            amr = run_isolates.amr
            serotype = run_isolates.serotype
            amr_status = amr |
                map { meta, amr_dir -> [ meta, "complete" ] }
            
        } else {
            amr = Channel.empty()
            mlst = Channel.empty()
            taxonomy_results = Channel.empty()
            sourmash_exclude = Channel.empty()
            serotype = Channel.empty()
            amr_status = reads |
                map { meta, reads, stats -> [ meta, "not-met" ] }
        }

        // Checkpoint 6 - AMR / isolates
        amr_checkpoint = amrCheckpoint(amr_status
        | mix( failed_samples )
        | unique() )

        bakta_version = baktaVersion()
        medaka_version = medakaVersion(bakta_version)
        mlst_version = mlstVersion(medaka_version)
        sourmash_version = sourmashVersion(mlst_version)
        mobsuite_version = mobsuiteVersion(sourmash_version)
        resfinder_version = resfinderVersion(mobsuite_version)
        seqsero_version = seqseroVersion(resfinder_version)
        software_versions = getVersions(seqsero_version)

        workflow_params = getParams()

        // Using meta.alias as join key, 
        // because directly joining meta objects introduces problems if some of the files are missing
        report_files_per_sample = reads 
            | filter {meta, reads, stats -> reads != null && meta != null }
            | map { meta, reads, stats_dir -> [meta.alias, meta, stats_dir ?: OPTIONAL_FILE] }
            | join(vcf_variant.map { meta, vcf -> [meta.alias, vcf] }, remainder: true)
            | join(vcf_stats.map { meta, stats -> [meta.alias, stats] }, remainder: true)
            | join(bakta.map { meta, gff3, gbff -> [meta.alias, [gff3, gbff]] }, remainder: true)
            | join(depth_stats.fwd.map { meta, bed -> [meta.alias, bed] }, remainder: true)
            | join(depth_stats.rev.map { meta, bed -> [meta.alias, bed] }, remainder: true)
            | join(depth_stats.all.map { meta, bed -> [meta.alias, bed] }, remainder: true)
            | join(flye_info.map { meta, stats -> [meta.alias, stats] }, remainder: true)
            | join(amr.map { meta, resfinder -> [meta.alias, resfinder] }, remainder: true)
            | join(mlst.map { meta, mlst_res -> [meta.alias, mlst_res] }, remainder: true)
            | join(serotype.map { meta, sero -> [meta.alias, sero] }, remainder: true)
            | join(taxonomy_results.map { meta, taxonomy -> [meta.alias, taxonomy] }, remainder: true)
            | join(mobsuite_results.map { meta, mob_res -> [meta.alias, mob_res] }, remainder: true)
            | combine(sourmash_exclude.ifEmpty([null]))
            | map { alias, meta, stats_dir, vcf_var, vcf_st, bakta_files, fwd, rev, all, flye, amr_res, mlst_res, sero, taxonomy, mob_res, excluded  ->
                def files = [stats_dir, vcf_var, vcf_st, flye, amr_res, mlst_res, taxonomy, sero, fwd, rev, all, mob_res, excluded]
                if (bakta_files) {
                    files.addAll(bakta_files)
                }
                [meta, files.findAll { it }]
                }

        sample_jsons = collect_results(report_files_per_sample, workflow_params)

        report_files_with_json = report_files_per_sample
            .join(sample_jsons, by: 0)
            .map { meta, files, json -> [meta, files + [json]] }

        run_model = createRunModel(
            sample_jsons.map { meta, json -> json }.collect(),
            metadata
        )

        fastq_stats = reads
        // replace `null` with path to optional file
        | map { [ it[0], it[1] ?: OPTIONAL_FILE, it[2] ?: OPTIONAL_FILE ] }
        | collectFastqIngressResultsInDir

        report = makeReport(
            software_versions,
            workflow_params,
            vcf_stats.map { meta, stats -> stats }.collect().ifEmpty(OPTIONAL_FILE),
            sample_ids.collect(),
            bakta_files.map { files -> params.run_bakta as Boolean && files.size() > 0 },
            fastq_stats.collect(),
            depth_stats.fwd.map{ meta, depths -> depths }.collect().ifEmpty(OPTIONAL_FILE),
            depth_stats.rev.map{ meta, depths -> depths }.collect().ifEmpty(OPTIONAL_FILE),
            depth_stats.all.map{ meta, depths -> depths }.collect().ifEmpty(OPTIONAL_FILE),
            run_model.collect(),
            client_fields)
        
        // Checkpoint 7 - report
        reporting_checkpoint = reportingCheckpoint(report)


        if (params.isolates) {
            perSampleReports = makePerSampleReports(
                software_versions,
                workflow_params,
                report_files_with_json
            )
            per_sample_report_status = perSampleReports 
                | map { meta, report -> [ meta, "complete" ] }
        } else {
            perSampleReports = Channel.empty()
            per_sample_report_status = reads |
                map { meta, reads, stats -> [ meta, "not-met" ] }
        }
        
        // Checkpoint 8 - per sample report
        per_sample_reporting_checkpoint = perSampleReportingCheckpoint(per_sample_report_status
        | mix( failed_samples )
        | unique() )

        accumulateCheckpoints.scan(
        ingress_checkpoint.mix(
            alignment_checkpoint,
            assembly_checkpoint,
            variant_checkpoint,
            plasmid_checkpoint,
            annotation_checkpoint,
            amr_checkpoint,
            reporting_checkpoint,
            per_sample_reporting_checkpoint
        ),
        metadata,
        definitions
    )

        all_out = vcf_stats.map{meta, stats -> stats}.concat(
            vcf_variant.map {meta, vcf -> vcf},
            consensus.map {meta, assembly -> assembly},
            alignments.map {meta, bam, bai -> [bam, bai]},
            report,
            perSampleReports.map {meta, report -> report},
            bakta.map{meta, gff3, gbff -> [gff3, gbff]},
            fastq_stats,
            amr.map {meta, resfinder -> resfinder},
            mlst.map {meta, mlst -> mlst},
            flye_info.map {meta, stats -> stats},
            workflow_params,
            software_versions,
            run_model,
            serotype.map { meta, sero -> sero },
            taxonomy_results.map {meta, taxonomy -> taxonomy},
            sourmash_exclude,
            mobsuite_results.map {meta, mob_res -> mob_res}
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

    data_root = "./data/"

    if (params.sourmash_db_exclude_list == null){
        params.remove('sourmash_db_exclude_list')
        params._sourmash_db_exclude_list = projectDir.resolve(data_root + "/sourmash_db_exclude_list.txt").toString()
    } else {
        params._sourmash_db_exclude_list = file(params.sourmash_db_exclude_list, type: "file", checkIfExists:true).toString()
        params.remove('sourmash_db_exclude_list')
    }

    if (params.sourmash_pointfinder_mapping == null){
        params.remove('sourmash_pointfinder_mapping')
        params._sourmash_pointfinder_mapping = projectDir.resolve(data_root + "/pointfinder_db_mapping.csv").toString()
    } else {
        params._sourmash_pointfinder_mapping = file(params.sourmash_pointfinder_mapping, type: "file", checkIfExists:true).toString()
        params.remove('sourmash_pointfinder_mapping')
    }

    // warn the user if overriding the basecall models found in the inputs
    if (params.override_basecaller_cfg) {
        log.warn \
            "Overriding basecall model with '${params.override_basecaller_cfg}'."
    }

    String fastcat_extra_args = params.min_read_length ? " -a $params.min_read_length" : ""

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
            "keep_unaligned":true,
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

