import groovy.json.JsonBuilder

process accumulateCheckpoints {
    label "wf_common"
    cpus 1
    memory "1 GB"
    input:
        path data
        val metadata
        path definitions
    output:
        path "checkpoints_${task.index}.json"
        val metadata
        path definitions
    publishDir params.out_dir, mode: 'copy', overwrite: true, pattern: "checkpoints_${task.index}.json", saveAs: { 'checkpoints.json' }
    script:
        // If the data list is lenth 1 then nextflow makes it not a list
        def data = data instanceof List ? data: [data]
        // Set-up the output file
        output = "checkpoints_${task.index}.json"
        // The run metadata (sample sheet)
        metaJson = new JsonBuilder(metadata).toPrettyString()
        // The checkpoint data
        checkpoint_data = data.getAt(0)
        // Our 1st checkpoint will not have a checkpoints file created
        // and the length of the data array wil be 1. Any subsequent checkpoint
        // will have the data followed by all of the previous checkpoint files
        // we need the last.
        if (data.size() > 1) {
            checkpoints_file = "--checkpoints_file ${data.getAt(-1)}"
        } else {
            checkpoints_file = ""
        }
    """
    echo '${metaJson}' > metadata.json
    accumulate_checkpoints.py ${output} \
        --output_definitions ${definitions} \
        --checkpoint_data ${checkpoint_data} \
        --metadata metadata.json \
        ${checkpoints_file}
    """
}

///////////////////////
// Workflow checkpoints
///////////////////////

process ingressCheckpoint {
  label "wf_common"
  cpus 1
  memory "1 GB"
  input:
    tuple val(meta), val(status)
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
     def files = status == 'complete' ? [ "sample-data" : "${meta.alias}/seqs.fastq.gz",
        "read-stats-per-file": "${meta.alias}/fastcat_stats/per-file-stats.tsv",
        "read-stats-per-read" : "${meta.alias}/fastcat_stats/per-read-stats.tsv.gz" ] : [ ]

    // make our checkpoint data
    def checkpoint_data = [[
            sample: "${meta.alias}",
            checkpoint_name: "reads",
            status: "${status}",
            files: files
        ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}

process alignmentCheckpoint {
  label "wf_common"
  cpus 1
  memory "1 GB"
  input:
    tuple val(meta), path(bam, stageAs: "bam/*"), path(bai, stageAs: "bai/*")
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
    def status = bam.fileName.name != "OPTIONAL_FILE" ? "complete" : "not-met"
    def files = status == 'complete' ? [ "alignment": "${bam.fileName.name}", "alignment-index": "${bai.fileName.name}" ] : [ ]

    // make our checkpoint data
    def checkpoint_data = [[
            sample: "${meta.alias}",
            checkpoint_name: "alignment",
            status: status,
            files: files
        ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}

process assemblyCheckpoint {
  label "wf_common"
  cpus 1
  memory "1 GB"
  input:
    tuple val(meta), path(fasta)
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
    def status = fasta.name != "OPTIONAL_FILE" ? "complete" : "not-met"
    def files = status == 'complete' ? [ "assembly": "${fasta}"] : [ ]

    // make our checkpoint data
    def checkpoint_data = [[
            sample: "${meta.alias}",
            checkpoint_name: "assembly",
            status: status,
            files: files
        ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}

process variantCheckpoint {
  label "wf_common"
  cpus 1
  memory "1 GB"
  input:
    tuple val(meta), val(status)
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
    def files = status == 'complete' ? [ "variants": "${meta.alias}.medaka.vcf.gz" ] : [ ]

    // make our checkpoint data
    def checkpoint_data = [[
            sample: "${meta.alias}",
            checkpoint_name: "variant-calling",
            status: "${status}",
            files: files
    ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}


process amrCheckpoint {
  label "wf_common"
  cpus 1
  memory "1 GB"
  input:
    tuple val(meta), val(status)
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
    def files = status == 'complete' ? [ "amr": "${meta.alias}_resfinder_results", "mlst": "${meta.alias}.mlst.json" ] : [ ]
    // make our checkpoint data
    def checkpoint_data = [[
            sample: "${meta.alias}",
            checkpoint_name: "amr",
            status: "${status}",
            files: files
        ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}


process annotationCheckpoint {
  label "wf_common"
  cpus 1
  memory "1 GB"
  input:
    tuple val(meta), val(status)
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
    def files = status == 'complete' ? [ "annotation": "${meta.alias}.prokka.gff" ] : [ ]
    // make our checkpoint data
    def checkpoint_data = [[        
            sample: "${meta.alias}",
            checkpoint_name: "annotation",
            status: "${status}",
            files: files
        ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}


process perSampleReportingCheckpoint {
  label "wf_common"
  cpus 1
  memory "1 GB"
  input:
    tuple val(meta), val(status)
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
    def files = status == 'complete' ? [ "per-sample-report": "${meta.alias}-isolate-report.html"  ] : [ ]
    def checkpoint_data = [[
            sample: "${meta.alias}",
            checkpoint_name: "per-sample-report",
            status: "${status}",
            files: files
        ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}


process reportingCheckpoint {
  label "wf_common"
  cpus 1
  memory "1 GB"
  input:
    path report
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
    String status = 'complete'

    // make our checkpoint data
    def checkpoint_data = [[
            sample: "",
            checkpoint_name: "reporting",
            status: "complete",
            files: [ 
              "workflow-report": "${report}",
              "workflow-results": "results.json"]
        ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}
