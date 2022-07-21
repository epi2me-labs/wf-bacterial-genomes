#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonBuilder

include { fastq_ingress } from './lib/fastqingress'
include { start_ping; end_ping } from './lib/ping'

process concatFastq {
    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple path(directory), val(meta)
    output:
        path "*reads.fastq"
        path "*stats*", emit: stats
        env SAMPLE_ID, emit: meta.sample_id
    shell:
    """
    # TODO: could do better here
    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} >  ${meta.sample_id}.reads.fastq
    SAMPLE_ID="${meta.sample_id}"
    """
}


process readStats {
    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_name), path(bam), path(bai)
    output:
        path "*readstats.txt", emit: stats
    """
    bamstats $bam > ${bam.simpleName}.readstats.txt
    if [[ \$(wc -l <${bam.simpleName}.readstats.txt) -le 1 ]]; then
        echo "No alignments of reads to reference sequence found."
        exit 1
    fi
    """
}


process coverStats {
    label "wfbacterialgenomes"
    cpus 2
    input:
        tuple val(sample_name), path(bam), path(bai)
    output:
        path "*fwd.regions.bed.gz", emit: fwd
        path "*rev.regions.bed.gz", emit: rev
        path "*total.regions.bed.gz", emit: all

    """
    mosdepth -n --fast-mode --by 200 --flag 16 -t $task.cpus ${sample_name}.fwd $bam
    mosdepth -n --fast-mode --by 200 --include-flag 16 -t $task.cpus ${sample_name}.rev $bam
    mosdepth -n --fast-mode --by 200 -t $task.cpus ${sample_name}.total $bam
    """
}


process deNovo {
    label "wfbacterialgenomes"
    cpus params.threads
    input:
        path reads
    output:
        path "*.fasta.gz"
    script:
        sample_name = reads.simpleName
        //mv output/00-assembly/draft_assembly.fasta ./${sample_name}.draft_assembly.fasta
    """
    flye --nano-raw $reads --genome-size $params.genome_size --out-dir output --threads $task.cpus
    mv output/assembly.fasta ./${sample_name}.draft_assembly.fasta
    bgzip ${sample_name}.draft_assembly.fasta
    """
}


process alignReads {
    label "wfbacterialgenomes"
    cpus params.threads
    input:
        tuple val(sample_name), path(reads), file(reference)
    output:
        tuple val("$sample_name"), path("*reads2ref.bam"), path("*reads2ref.bam.bai")
    script:
        sample_name = reads.simpleName
    """
    mini_align -i $reads -r $reference -p ${reads.simpleName}.reads2ref -t $task.cpus -m
    """
}


process splitRegions {
    // split the bam reference sequences into overlapping sub-regions

    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_name), path(bam), path(bai)
    output:
        stdout
    """
    #!/usr/bin/env python

    import itertools
    import medaka.common

    regions = itertools.chain.from_iterable(
        x.split($params.chunk_size, overlap=1000, fixed_size=False)
        for x in medaka.common.get_bam_regions("$bam"))
    region_list = []
    for reg in regions:
        print("$sample_name" + '.' + str(reg))
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
        tuple val(sample_name), val(reg), path(bam), path(bai)
    output:
        tuple val("$sample_name"), path("*consensus_probs.hdf")
    """
    medaka consensus $bam ${sample_name}.${task.index}.consensus_probs.hdf --threads 2 --model $params.medaka_model --region "$reg"
    """
}


process medakaVariant {

    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_name), path(consensus_hdf),  path(bam), path(bai), path(reference)
    output:
        path "*medaka.vcf.gz", emit: variants
        path '*.variants.stats', emit: variant_stats
    """
    medaka variant $reference *consensus_probs.hdf vanilla.vcf
    medaka tools annotate vanilla.vcf $reference $bam ${sample_name}.medaka.vcf
    bgzip -i ${sample_name}.medaka.vcf
    bcftools stats  ${sample_name}.medaka.vcf.gz > ${sample_name}.variants.stats
    """
}



process medakaConsensus {
    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_name), path(consensus_hdf),  path(bam), path(bai), path(reference)

    output:
        path "*medaka.fasta.gz"

    """
    medaka stitch --threads $task.cpus $consensus_hdf $reference ${sample_name}.medaka.fasta
    bgzip ${sample_name}.medaka.fasta
    """
}


process runProkka {
    // run prokka in a basic way on the consensus sequence
    label "prokka"
    cpus params.threads
    input:
        path consensus
    output:
        path "*prokka_results/*prokka.gbk"
    script:
        sample_name = "${consensus.simpleName}"

    """
    echo $sample_name
    gunzip -rf $consensus
    prokka ${params.prokka_opts} --outdir ${sample_name}.prokka_results --cpus $task.cpus --prefix ${sample_name}.prokka *medaka.fasta
    """
}


process getVersions {
    label "wfbacterialgenomes"
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
    return ch.map { it -> return tuple(it.split(/\./)[0], it.split(/\./)[1]) }
}

def nameIt(ch) {
    return ch.map { it -> return tuple(it.simpleName, it)}
}

def refTuple(ch) {
    return ch.map { it -> return tuple(it[0], it[2])}
}

// modular workflow
workflow calling_pipeline {
    take:
        reads
        reference
    main:


        reads = concatFastq(reads)
        if (!reference){
            ref  = deNovo(reads[0])
            println("No reference provided assuming De Novo.")
            named_refs = nameIt(ref)
            named_reads = nameIt(reads[0])
            refs_reads_groups = named_reads.join(named_refs)
            alignments = alignReads(refs_reads_groups)
            vcf_variant = Channel.empty()
            variants = Channel.empty()
        }
        else{
            references = channel.fromPath(params.reference)
            named_reads = nameIt(reads[0])
            ref_reads_groups = named_reads.combine(references)
            alignments = alignReads(ref_reads_groups)
            named_refs = refTuple(ref_reads_groups)

        }

        sample_ids = reads.sample_id.collectFile(name: 'sample_ids.csv', newLine: true)
        read_stats = readStats(alignments)
        depth_stats = coverStats(alignments)
        regions = splitRegions(alignments).splitText()
        named_regions = groupIt(regions)
        regions_bams = named_regions.combine(alignments, by: [0])
        hdfs = medakaNetwork(regions_bams)
        hdfs_grouped = hdfs.groupTuple().combine(alignments, by: [0]).join(named_refs)
        consensus = medakaConsensus(hdfs_grouped)
        if (reference){
            variant = medakaVariant(hdfs_grouped)
            variants = variant.variant_stats
            vcf_variant = variant.variants
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
    emit:
        vcf_variant
        variants
        consensus
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

    if (!params.fastq) {
        helpMessage()
        println("")
        println("`--fastq` is required and --reference is required when performing variant calling")
        exit 1
    }
    samples = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet, params.sanitize_fastq)
    reference = params.reference
    results = calling_pipeline(samples, reference)
    output(results.consensus.concat(results.report, results.prokka, results.variants))
    end_ping(calling_pipeline.out.telemetry)

}
