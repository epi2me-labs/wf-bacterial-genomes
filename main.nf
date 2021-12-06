#!/usr/bin/env extflow

nextflow.enable.dsl = 2
import groovy.json.JsonBuilder

include { fastq_ingress } from './lib/fastqingress'
include { start_ping; end_ping } from './lib/ping'

process concatFastq {
    label "containerCPU"
    cpus 1
    input:
        tuple path(directory), val(sample_id), val(type)
    output:
        file "reads.fastq"

    shell:
    """
    # TODO: could do better here
    fastcat -s ${sample_id} -r ${sample_id}.stats -x ${directory} >  reads.fastq
    """
}


process readStats {
    label "containerCPU"
    cpus 1
    input:
        tuple file("alignments.bam"), file("alignments.bam.bai")
    output:
        path "readstats.txt", emit: stats
        path "mapsummary.txt", emit: summary
    """
    stats_from_bam -s mapsummary.txt -o readstats.txt alignments.bam
    if [[ \$(wc -l <readstats.txt) -le 1 ]]; then
        echo "No alignments of reads to reference sequence found."
        exit 1
    fi
    """
}


process coverStats {
    label "containerCPU"
    cpus 2
    input:
        tuple file("alignments.bam"), file("alignments.bam.bai")
    output:
        file "depth.txt"
    """
    coverage_from_bam --one_file depth.txt --stride 100 alignments.bam
    """
}


process alignReads {
    label "containerCPU"
    cpus params.threads
    input:
        file reads
        file reference
    output:
        tuple file("reads2ref.bam"), file("reads2ref.bam.bai")

    """
    mini_align -i $reads -r $reference -p reads2ref -t $task.cpus -m
    """
}


process splitRegions {
    // split the bam reference sequences into overlapping sub-regions

    label "containerCPU"
    cpus 1
    input:
        tuple file(bam), file(bai)
    output:
        stdout
    """
    #!/usr/bin/env python

    import itertools
    import medaka.common

    regions = itertools.chain.from_iterable(
        x.split($params.chunk_size, overlap=1000, fixed_size=False)
        for x in medaka.common.get_bam_regions("$bam"))
    for reg in regions:
        print(reg)
    """
}


// TODO: in a single GPU environment it would be better just
//       to use a single process for the whole bam file. Need
//       to read up on conditional channels

process medakaNetwork {
    // run medaka consensus for each region

    label "containerCPU"
    cpus 2
    input:
        tuple file(bam), file(bai)
        each reg
    output:
        file "consensus_probs.hdf"
    """
    medaka consensus $bam consensus_probs.hdf --threads 2 --model $params.medaka_model --region "$reg"
    """
}


process medakaVariant {
    label "containerCPU"
    cpus 1
    input:
        file "consensus_probs_*.hdf"
        tuple file(bam), file(bai)
        file reference
    output:
        tuple path("medaka.vcf.gz"), path("medaka.vcf.gz.gzi")
    """
    medaka variant $reference consensus_probs_*.hdf vanilla.vcf
    medaka tools annotate vanilla.vcf $reference $bam medaka.vcf
    bgzip -i medaka.vcf
    """
}


process medakaConsensus {
    label "containerCPU"
    cpus 1
    input:
        file "consensus_probs_*.hdf"
        file reference
    output:
        file "medaka.fasta.gz"

    """
    medaka stitch --threads $task.cpus consensus_probs_*.hdf $reference medaka.fasta
    bgzip medaka.fasta
    """
}


process runProkka {
    // run prokka in a basic way on the consensus sequence
    label "prokka"
    cpus 1
    input:
        file "consensus.fasta"
    output:
        file "prokka_results"
        
    """
    prokka ${params.prokka_opts} --outdir prokka_results --prefix prokka consensus.fasta

    """
}


process getVersions {
    label "containerCPU"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    medaka --version | sed 's/ /,/' >> versions.txt
    python -c "import pomoxis; print(f'pomoxis,{pomoxis}.__version__')" >> versions.txt
    python -c "import tensorflow; print(f'tensorflow,{tensorflow}.__version__')" >> versions.txt
    """
}


process getParams {
    label "containerCPU"
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
    label "containerCPU"
    cpus 1
    input:
        file "depth.txt"
        file "read_summary.txt"
        file "align_summary.txt"
        tuple path("medaka.vcf.gz"), path("medaka.vcf.gz.gzi")
        path "versions/*"
        path "params.json"

    output:
        file "wf-hap-snps-*.html"
    script:
        report_name = "wf-hap-snps-" + params.report_name + '.html'
    """
    bcftools stats medaka.vcf.gz > variants.stats
    report.py depth.txt read_summary.txt align_summary.txt variants.stats $report_name \
        --versions versions \
        --params params.json

    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "containerCPU"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


// modular workflow
workflow calling_pipeline {
    take:
        reads
        reference
    main:
        reads = concatFastq(reads)
        alignments = alignReads(reads, reference)
        read_stats = readStats(alignments)
        depth_stats = coverStats(alignments)
        regions = splitRegions(alignments).splitText()
        hdfs = medakaNetwork(alignments, regions).collect()
        consensus = medakaConsensus(hdfs, reference)
        variants = medakaVariant(hdfs, alignments, reference)
        if (params.run_prokka) {
            prokka = runProkka(consensus)
        } else {
            prokka = Channel.empty()
        } 
        software_versions = getVersions()
        workflow_params = getParams()
        report = makeReport(depth_stats, read_stats.stats, read_stats.summary, variants,
                            software_versions.collect(), workflow_params)
        telemetry = workflow_params
    emit:
        consensus
        variants
        report
        prokka
        telemetry
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    start_ping()
    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq or !params.reference) {
        helpMessage()
        println("")
        println("`--fastq` and `--reference` are required")
        exit 1
    }
    samples = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet, params.sanitize_fastq)
    reference = channel.fromPath(params.reference)
    results = calling_pipeline(samples, reference)
    output(results.consensus.concat(results.variants, results.report, results.prokka))
    end_ping(calling_pipeline.out.telemetry)

}
