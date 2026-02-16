Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | wf-bacterial-genomes-report.html | Report for all samples. | aggregated |
| Workflow results | results.json | Structured workflow results for internal/onward use. | aggregated |
| Workflow checkpoints | checkpoints.json | Structured workflow checkpoints for internal/onward use. | aggregated |
| Concatenated sequence data | {{ alias }}/seqs.fastq.gz | Per sample reads concatenated in to one fastq file. | per-sample |
| Per file read stats | {{ alias }}/fastcat_stats/per-file-stats.tsv | A TSV with per file read stats. | per-sample |
| Per read stats | {{ alias }}/fastcat_stats/per-read-stats.tsv.gz | A TSV with per read stats. | per-sample |
| Flye assembly stats | {{ alias }}.flye_stats.tsv | Assembly statistics generated for Flye assemblies. | per-sample |
| Alignment | {{ alias }}.bam | Aligned reads for the sample in BAM format. | per-sample |
| Alignment index | {{ alias }}.bam.bai | An index file for the alignment in BAI format. | per-sample |
| Draft assembly FASTA file | {{ alias }}.medaka.fasta.gz | Consensus file generated from either de novo assembly or reference variant calling pipeline, with contigs oriented to start at replication initiation point. | per-sample |
| Plasmid contig identification results | {{ alias }}_mob_results | Complete MOB-suite output directory containing all plasmid identification and characterisation results | per-sample |
| Variants VCF file | {{ alias }}.medaka.vcf.gz | VCF file of variants detected against the provided reference (reference mode only). | per-sample |
| Variants summary | {{ alias }}.variants.stats | TSV file of summary statistics for variants in sample (reference mode only). | per-sample |
| Annotations in GFF3 format | {{ alias }}.bakta.gff3 | Annotations of regions of interest in assembly in GFF3 format. | per-sample |
| Annotations in GBFF format | {{ alias }}.bakta.gbff | Annotations of regions of interest in assembly in GBFF format. | per-sample |
| Sequence typing results | {{ alias }}.mlst.json | Sequence typing results in JSON format (isolates mode only). | per-sample |
| Sourmash taxonomic identification results with lineage annotations. | {{ alias }}_sourmash_taxonomy.csv | Sourmash gather hits with taxonomic lineage annotations (isolates mode only). | per-sample |
| Sourmash excluded assemblies | sourmash_picklist_excluded.txt | List of assembly identifiers (GCA_/GCF_) excluded from default Sourmash species identification database (isolates mode only). | aggregated |
| AMR calling results | {{ alias }}_resfinder_results | Complete AMR calling output directory containing detailed ResFinder, PointFinder, and DesinFinder outputs (isolates mode only). | per-sample |
| Isolate per sample report | {{ alias }}-isolates-report.html | Per sample report (isolates mode only). | per-sample |
