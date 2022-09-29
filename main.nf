#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonBuilder

include { fastq_ingress } from './lib/fastqingress'


process concatFastq {
    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple path(directory), val(meta)
    output:
        tuple val(meta.sample_id), path("*reads.fastq"), emit: read
        path "*stats*", emit: stats
        env SAMPLE_ID, emit: sample_id
    shell:
    """
    # TODO: could do better here
    fastcat -s "${meta.sample_id}" -r "${meta.sample_id}.stats" -x "${directory}" >  "${meta.sample_id}.reads.fastq"
    SAMPLE_ID="${meta.sample_id}"
    """
}


process readStats {
    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_id), path(bam), path(bai)
    output:
        path "*readstats.txt", emit: stats
    """
    bamstats "${bam}" > "${sample_id}.readstats.txt"
    if [[ \$(wc -l <"${sample_id}.readstats.txt") -le 1 ]]; then
        echo "No alignments of reads to reference sequence found."
        exit 1
    fi
    """
}


process coverStats {
    label "wfbacterialgenomes"
    cpus 2
    input:
        tuple val(sample_id), path(bam), path(bai)
    output:
        path "*fwd.regions.bed.gz", emit: fwd
        path "*rev.regions.bed.gz", emit: rev
        path "*total.regions.bed.gz", emit: all

    """
    mosdepth -n --fast-mode --by 200 --flag 16 -t $task.cpus "${sample_id}.fwd" "${bam}"
    mosdepth -n --fast-mode --by 200 --include-flag 16 -t $task.cpus "${sample_id}.rev" "${bam}"
    mosdepth -n --fast-mode --by 200 -t $task.cpus "${sample_id}.total" "${bam}"
    """
}


process deNovo {
    label "wfbacterialgenomes"
    cpus params.threads
    input:
        tuple val(sample_id), path(reads)
    output:
        tuple val(sample_id), path("${sample_id}.draft_assembly.fasta.gz")
    """
    flye --nano-raw "${reads}" --genome-size "${params.genome_size}" --out-dir output --threads "${task.cpus}"
    mv output/assembly.fasta "./${sample_id}.draft_assembly.fasta"
    bgzip "${sample_id}.draft_assembly.fasta"
    """
}


process alignReads {
    label "wfbacterialgenomes"
    cpus params.threads
    input:
        tuple val(sample_id), path(reads), file(reference)
    output:
        tuple val(sample_id), path("*reads2ref.bam"), path("*reads2ref.bam.bai")
    """
    mini_align -i "${reads}" -r "${reference}" -p "${sample_id}.reads2ref" -t $task.cpus -m
    """
}


process splitRegions {
    // split the bam reference sequences into overlapping sub-regions

    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_id), path(bam), path(bai)
    output:
        stdout
    """
    #!/usr/bin/env python

    import itertools
    import medaka.common

    regions = itertools.chain.from_iterable(
        x.split(${params.chunk_size}, overlap=1000, fixed_size=False)
        for x in medaka.common.get_bam_regions("${bam}"))
    region_list = []
    for reg in regions:
        print("${sample_id}" + '&split!' + str(reg))
    """
}


// TODO: in a single GPU environment it would be better just
//       to use a single process for the whole bam file. Need
//       to read up on conditional channels

process medakaNetwork {
    // run medaka consensus for each region

    label "wfbacterialgenomes"
    cpus 2
    input:
        tuple val(sample_id), val(reg), path(bam), path(bai)
    output:
        tuple val(sample_id), path("*consensus_probs.hdf")
    """
    medaka consensus "${bam}" "${sample_id}.${task.index}.consensus_probs.hdf" \
        --threads 2 --model "${params.medaka_model}" --region "${reg}"
    """
}


process medakaVariant {

    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_id), path(consensus_hdf),  path(bam), path(bai), path(reference)
    output:
        path "${sample_id}.medaka.vcf.gz", emit: variants
        path "${sample_id}.variants.stats", emit: variant_stats
    """
    medaka variant "${reference}" *consensus_probs.hdf vanilla.vcf
    medaka tools annotate vanilla.vcf "${reference}" "${bam}" "${sample_id}.medaka.vcf"
    bgzip -i "${sample_id}.medaka.vcf"
    bcftools stats  "${sample_id}.medaka.vcf.gz" > "${sample_id}.variants.stats"
    """
}



process medakaConsensus {
    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_id), path(consensus_hdf),  path(bam), path(bai), path(reference)
    output:
        tuple val(sample_id), path("${sample_id}.medaka.fasta.gz")

    """
    medaka stitch --threads $task.cpus "${consensus_hdf}" "${reference}" "${sample_id}.medaka.fasta"
    bgzip "${sample_id}.medaka.fasta"
    """
}


process runProkka {
    // run prokka in a basic way on the consensus sequence
    label "prokka"
    cpus params.threads
    input:
        tuple val(sample_id), path(consensus)
    output:
        path "*prokka_results/*prokka.gbk"
    script:
        def prokka_opts = "${params.prokka_opts}" == null ? "${params.prokka_opts}" : ""
    """
    echo $sample_id
    gunzip -rf $consensus
    prokka $prokka_opts --outdir "${sample_id}.prokka_results" \
        --cpus $task.cpus --prefix "${sample_id}.prokka" *medaka.fasta
    """
}


process getVersions {
    label "wfbacterialgenomes"
    cpus 1
    output:
        path "versions.txt"
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    medaka --version | sed 's/ /,/' >> versions.txt
    python -c "import pomoxis; print(f'pomoxis,{pomoxis}.__version__')" >> versions.txt
    python -c "import tensorflow; print(f'tensorflow,{tensorflow}.__version__')" >> versions.txt
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
        file "params.json"
        path "variants/*"
        path sample_ids
        path "prokka/*"
        path "stats/*"
        path "fwd/*"
        path "rev/*"
        path "total_depth/*"
    output:
        path "wf-bacterial-genomes-*.html"
    script:
        report_name = "wf-bacterial-genomes-" + params.report_name + '.html'
        prokka = params.run_prokka as Boolean ? "--prokka prokka/*" : ""
    """
    report.py --bcf_stats variants/* \
    $prokka \
    --versions versions \
    --params params.json \
    --sample_ids $sample_ids \
    --output $report_name \
    --stats stats/*

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
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


def groupIt(ch) {
    return ch.map { it -> return tuple(it.split(/&split!/)[0], it.split(/&split!/)[1]) }
}


// modular workflow
workflow calling_pipeline {
    take:
        reads
        reference
    main:
        reads = concatFastq(reads)
        if (!reference){
            println("No reference provided creating de-novo assemblies.")
            named_refs = deNovo(reads.read)
            read_ref_groups = reads.read.join(named_refs)
        } else {
            references = channel.fromPath(params.reference)
            read_ref_groups = reads.read.combine(references)
            named_refs = read_ref_groups.map { it -> [it[0], it[2]] }
        }
        read_ref_groups.view()
        alignments = alignReads(read_ref_groups)

        sample_ids = reads.sample_id.collectFile(name: 'sample_ids.csv', newLine: true)
        read_stats = readStats(alignments)
        depth_stats = coverStats(alignments)
        regions = splitRegions(alignments).splitText()
        named_regions = groupIt(regions)
        
        regions_bams = named_regions.combine(alignments, by: [0])
        hdfs = medakaNetwork(regions_bams)
        hdfs_grouped = hdfs.groupTuple().combine(alignments, by: [0]).join(named_refs)
        consensus = medakaConsensus(hdfs_grouped)

        // call variants
        if (reference){
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

        software_versions = getVersions()
        workflow_params = getParams()

        report = makeReport(
                            software_versions.collect(),
                            workflow_params,
                            variants.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
                            sample_ids,
                            prokka.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
                            reads.stats.collect(),
                            depth_stats.fwd.collect(),
                            depth_stats.rev.collect(),
                            depth_stats.all.collect())
        telemetry = workflow_params
        all_out = variants.concat(vcf_variant,
                      consensus.map {it -> it[1]},
                      report,
                      prokka)
        
    emit:
        all_out
        telemetry
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    if (params.disable_ping == false) {
        try { 
            Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
        } catch(RuntimeException e1) {
        }
    }
    
    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq) {
        helpMessage()
        println("")
        println("`--fastq` is required and --reference is required when performing variant calling")
        exit 1
    }

    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "sanitize": params.sanitize_fastq,
        "output":params.out_dir])

    reference = params.reference
    results = calling_pipeline(samples, reference)
    output(results.all_out)
    

}

if (params.disable_ping == false) {
    workflow.onComplete {
        try{
            Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
        }catch(RuntimeException e1) {
        }
    }
    
    workflow.onError {
        try{
            Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
        }catch(RuntimeException e1) {
        }
    }
}
