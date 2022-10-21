#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonBuilder

include { fastq_ingress } from './lib/fastqingress'


process concatFastq {
    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple path("input"), val(meta)
    output:
        tuple val(meta.sample_id), path("${meta.sample_id}.reads.fastq.gz"), emit: read
        path "*stats*", emit: stats
    shell:
    """
    # TODO: could do better here
    fastcat -s "${meta.sample_id}" -r "${meta.sample_id}.stats" -x input | bgzip > "${meta.sample_id}.reads.fastq.gz"
    SAMPLE_ID="${meta.sample_id}"
    """
}


process readStats {
    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_id), path("align.bam"), path("align.bam.bai")
    output:
        path "*readstats.txt", emit: stats
    """
    bamstats align.bam > "${sample_id}.readstats.txt"
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
        tuple val(sample_id), path("align.bam"), path("align.bam.bai")
    output:
        path "*fwd.regions.bed.gz", emit: fwd
        path "*rev.regions.bed.gz", emit: rev
        path "*total.regions.bed.gz", emit: all

    """
    mosdepth -n --fast-mode --by 200 --flag 16 -t $task.cpus "${sample_id}.fwd" align.bam
    mosdepth -n --fast-mode --by 200 --include-flag 16 -t $task.cpus "${sample_id}.rev" align.bam
    mosdepth -n --fast-mode --by 200 -t $task.cpus "${sample_id}.total" align.bam
    """
}


process deNovo {
    label "wfbacterialgenomes"
    cpus params.threads
    input:
        tuple val(sample_id), path("reads.fastq.gz")
    output:
        tuple val(sample_id), path("${sample_id}.draft_assembly.fasta.gz"), path("${sample_id}_flye_stats.tsv")
        
    """
    flye --nano-raw reads.fastq.gz --genome-size "${params.genome_size}" --out-dir output --threads "${task.cpus}"
    mv output/assembly.fasta "./${sample_id}.draft_assembly.fasta"
    mv output/assembly_info.txt "./${sample_id}_flye_stats.tsv"
    bgzip "${sample_id}.draft_assembly.fasta"
    """
}


process alignReads {
    label "wfbacterialgenomes"
    cpus params.threads
    input:
        tuple val(sample_id), path("reads.fastq.gz"), path("ref.fasta.gz")
    output:
        tuple val(sample_id), path("*reads2ref.bam"), path("*reads2ref.bam.bai")
    """
    mini_align -i reads.fastq.gz -r ref.fasta.gz -p "${sample_id}.reads2ref" -t $task.cpus -m
    """
}


process splitRegions {
    // split the bam reference sequences into overlapping sub-regions

    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_id), path("align.bam"), path("align.bam.bai")
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
        tuple val(sample_id), val(reg), path("align.bam"), path("align.bam.bai")
    output:
        tuple val(sample_id), path("*consensus_probs.hdf")
    """
    medaka consensus align.bam "${sample_id}.consensus_probs.hdf" \
        --threads 2 --model "${params.medaka_model}" --region "${reg}"
    """
}


process medakaVariant {
    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_id), path("consensus_probs*.hdf"),  path("align.bam"), path("align.bam.bai"), path("ref.fasta.gz")
    output:
        path "${sample_id}.medaka.vcf.gz", emit: variants
        path "${sample_id}.variants.stats", emit: variant_stats
    // note: extension on ref.fasta.gz might not be accurate but shouldn't (?) cause issues.
    //       Also the first step may create an index if not already existing so the alternative
    //       reference.* will break 
    """
    medaka variant ref.fasta.gz consensus_probs*.hdf vanilla.vcf
    medaka tools annotate vanilla.vcf ref.fasta.gz align.bam "${sample_id}.medaka.vcf"
    bgzip -i "${sample_id}.medaka.vcf"
    bcftools stats  "${sample_id}.medaka.vcf.gz" > "${sample_id}.variants.stats"
    """
}

process assemblyStats {
    label "wfbacterialgenomes"
    input:
         path(sample_assembly_gz)

    output:
        path("quast_output/transposed_report.tsv")

    """
    quast -o quast_output -t $task.cpus ${sample_assembly_gz}

    """
}



process medakaConsensus {
    label "wfbacterialgenomes"
    cpus 1
    input:
        tuple val(sample_id), path("consensus_probs*.hdf"),  path("align.bam"), path("align.bam.bai"), path("reference*")
    output:
        tuple val(sample_id), path("${sample_id}.medaka.fasta.gz")

    """
    medaka stitch --threads $task.cpus consensus_probs*.hdf reference* "${sample_id}.medaka.fasta"
    bgzip "${sample_id}.medaka.fasta"
    """
}


process runProkka {
    // run prokka in a basic way on the consensus sequence
    label "prokka"
    cpus params.threads
    input:
        tuple val(sample_id), path("consensus.fasta.gz")
    output:
        path "*prokka_results/*prokka.gbk"
    script:
        def prokka_opts = "${params.prokka_opts}" == null ? "${params.prokka_opts}" : ""
    """
    echo $sample_id
    gunzip -rf consensus.fasta.gz
    prokka $prokka_opts --outdir "${sample_id}.prokka_results" \
        --cpus $task.cpus --prefix "${sample_id}.prokka" *consensus.fasta
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
        path "params.json"
        path "variants/*"
        val sample_ids
        path "prokka/*"
        path "stats/*"
        path "fwd/*"
        path "rev/*"
        path "total_depth/*"
        path "assembly_QC.txt"
        path "flye_stats/*"
    output:
        path "wf-bacterial-genomes-*.html"
    script:
        report_name = "wf-bacterial-genomes-" + params.report_name + '.html'
        prokka = params.run_prokka as Boolean ? "--prokka" : ""
        denovo = params.reference == null ? "--denovo" : ""
        samples = sample_ids.join(" ")
    // NOTE: the script assumes the various subdirectories
    """
    report.py \
    $prokka $denovo \
    --versions versions \
    --params params.json \
    --output $report_name \
    --sample_ids $samples \
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


// modular workflow
workflow calling_pipeline {
    take:
        reads
        reference
    main:
        reads = concatFastq(reads)
        sample_ids = reads.read.map { it -> it[0] }
        if (!reference){
            println("No reference provided creating de-novo assemblies.")
            denovo_assem = deNovo(reads.read)
            named_refs = denovo_assem.map { it -> [it[0], it[1]] }
            read_ref_groups = reads.read.join(named_refs)
            

        } else {
            references = channel.fromPath(params.reference)
            read_ref_groups = reads.read.combine(references)
            named_refs = read_ref_groups.map { it -> [it[0], it[2]] }
        }
        alignments = alignReads(read_ref_groups)

        read_stats = readStats(alignments)
        depth_stats = coverStats(alignments)
        regions = splitRegions(alignments).splitText()
        named_regions = regions.map {
            it -> return tuple(it.split(/&split!/)[0], it.split(/&split!/)[1])
        }
        
        regions_bams = named_regions.combine(alignments, by: [0])
        hdfs = medakaNetwork(regions_bams)
        hdfs_grouped = hdfs.groupTuple().combine(alignments, by: [0]).join(named_refs)
        consensus = medakaConsensus(hdfs_grouped)

        // post polishing, do assembly specific things
        if (!reference){
             println("No reference provided, analysing assemblies.")
             assem_stats = assemblyStats(consensus.collect({it -> it[1]}))
             flye_info = denovo_assem.map { it -> it[2] }
        } else {
             assem_stats = Channel.empty()
             flye_info = Channel.empty()
        }

             

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
            variants.collect().ifEmpty(file("${projectDir}/data/OPTIONAL_FILE")),
            sample_ids.collect(),
            prokka.collect().ifEmpty(file("${projectDir}/data/OPTIONAL_FILE")),
            reads.stats.collect(),
            depth_stats.fwd.collect(),
            depth_stats.rev.collect(),
            depth_stats.all.collect(),
            assem_stats.ifEmpty(file("${projectDir}/data/OPTIONAL_FILE")),
            flye_info.collect().ifEmpty(file("${projectDir}/data/OPTIONAL_FILE")))
        telemetry = workflow_params
        all_out = variants.concat(
            vcf_variant,
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
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }
    
    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet])

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
