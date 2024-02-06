Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | wf-bacterial-genomes-report.html | Report for all samples | aggregated |
| Workflow results | results.json | Structured workflow results for internal/onward use. | aggregated |
| workflow checkpoints | checkpoints.json | Structured workflow checkpoints for internal/onward use. | aggregated |
| Concatenated sequence data | {{ alias }}/seqs.fastq.gz | Per sample reads concatenated in to one fastq file. | per-sample |
| Per file read stats | {{ alias }}/fastcat_stats/per-file-stats.tsv | A TSV with per file read stats, per sample. | per-sample |
| Per read stats | {{ alias }}/fastcat_stats/per-read-stats.tsv.gz | A TSV with per read stats, per sample | per-sample |
| alignment | {{ alias }}.bam | Aligned reads for the sample in BAM format. | per-sample |
| alignment index | {{ alias }}.bam.bai | An index file for the alignment in BAI format. | per-sample |
| Draft assembly FASTA file | {{ alias }}.medaka.fasta.gz | Consensus file generated from either de-novo assembly or reference variant calling pipeline. | per-sample |
| Variants VCF file | {{ alias }}.medaka.vcf.gz | VCF file of variants detected against the provided reference (Reference mode only). | per-sample |
| Variants summary | {{ alias }}.variants.stats | TSV file of summary statistics for variants in sample (Reference mode only). | per-sample |
| Annotations files | {{ alias }}.prokka.gff | Annotations of regions of interest in assembly in GBK and GFF format. | per-sample |
| Sequence typing results | {{ alias }}.mlst.json | Sequence typing results in JSON format (isolates mode only). | per-sample |
| AMR calling results | {{ alias }}_resfinder_results | Resfinder results for AMR calling (isolates mode only). | per-sample |
| isolates per sample report | {{ alias }}-isolates-report.html | Per sample report isolates mode | per-sample |
