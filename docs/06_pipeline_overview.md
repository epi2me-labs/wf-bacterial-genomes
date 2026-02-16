<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->

### 1. Preprocessing concatenates input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read statistics, such as read length, count and quality statistics.

### 2a. De novo assembly

#### i. Assembly

[Flye](https://github.com/fenderglass/Flye) is used to generate draft assembly from the FASTQ reads. By default, Flye runs with `--nano-hq` parameter. Additional configuration can be specified using `--flye_opts` parameter.

Flye produces a draft assembly along with a statistics file containing information on contig lengths, coverage, multiplicity, and circularity. For de novo assemblies, contigs are assigned random names, whereas in reference-based assemblies, contig names are inherited from the reference sequence.

#### ii. Contig reorientation 

Following de novo assembly, circular contigs are processed with [dnaapler](https://github.com/gbouras13/dnaapler), which reorients contigs to a common start position at the origin of replication: _dnaA_ for chromosomes and _repA_ for plasmids.
Dnaapler relies on MMseqs2-based search against a database consisting of origin of replication sequences.
Information on the database content can be found in the [dnaapler documentation](https://github.com/gbouras13/dnaapler).

#### iii. Alignment

Following contig reorientation, reads are aligned against the assembly using [mini_align](https://github.com/nanoporetech/pomoxis/).
Subsequently, [bamstats](https://github.com/epi2me-labs/fastcat?tab=readme-ov-file#bamstats) and [mosdepth](https://github.com/brentp/mosdepth) are used to calculate mapped read and genome coverage statistics.

#### iv. Polishing

The draft assembly from Flye is then polished using [medaka](https://github.com/nanoporetech/medaka).
This step will attempt to correct any errors that were introduced during the de novo assembly process and generate a FASTA output file.

The workflow selects the appropriate [medaka models](https://github.com/nanoporetech/medaka#models) based on the basecaller configuration that was used to process the signal data.
Per default, the workflow will attempt to determine the basecaller model from the input data.
When this fails (or when you wish to override the automatic selection), it can be provided with `--override_basecaller_cfg`.

### 2b. Variant calling mode

#### i. Align reads

Reads are aligned against the provided reference with [mini_align](https://github.com/nanoporetech/pomoxis/).

#### ii. Call variants

After alignment, haploid variants are called with [medaka](https://github.com/nanoporetech/medaka) (see the Polishing section above for details on Medaka model selection).

#### iii. Use the variants to generate a consensus

The variants passing the depth filter are then incorporated in the reference to create the consensus sequence. Variant stats are also created at this point.

### 3. Plasmid identification

Plasmid identification and characterization is performed using [MOB-suite](https://github.com/phac-nml/mob-suite).
This workflow implements `mob-recon`, which is used to perform plasmid reconstruction and typing.
Contigs are characterized as chromosomal or plasmid.
Plasmid contigs sharing similar replication and mobilization characteristics are grouped into MOB-clusters, representing distinct plasmids that may consist of multiple fragments.

The summary report includes a section listing plasmid contigs and identified characteristics, such as replicon, relaxase, MPF, and _oriT_ typing information, as well as mobility prediction (conjugative, mobilizable, or non-mobilizable), and predicted/observed plasmid host range.
The output directory contains a `{alias}_mob_results` subdirectory, which includes a `mobtyper-results.txt` file with additional characteristics for each contig.
If plasmid contigs are identified, this directory will also contain `.fasta` files grouping plasmid contigs based on MOB-cluster identification, and an additional `mobtyper-results.txt` file with extended typing information.

Plasmid identification can be disabled using the `--run_plasmid_id` parameter.
Additional [mob-recon options](https://github.com/phac-nml/mob-suite) can be specified using the `--plasmid_id_opts` parameter.

### 4. Annotations

Consensus genome assemblies can be annotated using [Bakta](https://github.com/oschwengers/bakta). 
By default, Bakta will run with [light database](https://github.com/oschwengers/bakta?tab=readme-ov-file#database), which will be automatically downloaded and installed when the workflow is executed. To optionally run Bakta with the [full database](https://zenodo.org/records/14916843) specify `--bakta_db_type=full`. This will trigger the installation of the full version of the database. As of version 6.0, the compressed size of the full database is approximately 32 GB, so the download may take a significant amount of time.
Alternatively, users can use following steps to manually download and install Bakta database, and specify the custom database location using the `--bakta_db` parameter.

#### Downloading the full database (via Docker)
A compatible database version can be downloaded using the following command:

```bash
docker run -v /path/to/desired-db-path:/db --entrypoint /bin/bash ontresearch/bakta:latest
 -c "bakta_db download --output /db --type full"
```

#### Manual download

Alternatively, you can manually download the database:
```bash
wget https://zenodo.org/record/14916843/files/db.tar.xz
```

Then install it using Docker:
```bash
docker run -v /path/to/db:/db --entrypoint /bin/bash ontresearch/bakta:latest -c "bakta_db install -i db.tar.xz"
```

#### Updating AMRFinderPlus database
Bakta relies on the AMRFinderPlus database, which may need to be updated manually (even after a fresh installation of the Bakta database). To update it:

```bash
docker run -v /path/to/db:/db --entrypoint /bin/bash ontresearch/bakta:latest -c "amrfinder_update --force_update --database db/amrfinderplus-db/"
```

#### Using a local database in the workflow
To use a locally installed Bakta database, pass the path using `--bakta_db`. For example:

```bash
nextflow run epi2me-labs/wf-bacterial-genomes \
	--fastq '/path/to/fastq' \
	--bakta_db '/path/to/db'
```

By default, the workflow disables circular genome plot generation by Bakta (`--skip-plot`), but all other settings follow Bakta’s defaults. Users can customize Bakta's behavior by supplying additional arguments via `--bakta_opts`.

### 5. Isolates mode (optional)

#### i. Multi-locus sequence typing (MLST)

MLST is a common technique used to help characterise bacterial isolates by using allelic variation from internal DNA fragments of 6-7 housekeeping genes. Typing schemes for specific species and genera are found on [PubMLST](https://pubmlst.org/) and are pre-loaded into this workflow. [MLST](https://github.com/tseemann/mlst) is used to scan assembly and automatically infer the correct typing scheme, and subsequently identify the allele variant found.

MLST results determine the species assignment shown in the final report. This information also determines whether the sample will undergo Salmonella typing or point mutation–based antimicrobial resistance prediction.

#### ii. Species identification

[Sourmash](https://sourmash.readthedocs.io/en/latest/index.html) is used to perform taxonomic identification of assemblies.
Sourmash performs a k-mer based sequence comparison against genomes in a reference database collection.
The analysis uses MinHash sketching to rapidly identify the most similar reference genomes and provide species-level classification with confidence metrics.
Sourmash results provide taxonomic assignment from domain to species level, along with Average Nucleotide Identity (ANI) values and coverage statistics, to assess the quality of the match.
The analysis identifies the primary species match as well as secondary matches that may represent genomic variation within the species.
Species identification complements MLST analysis by providing broader taxonomic context and can help resolve cases where MLST schemes may not be available for the organism.
The isolate HTML report includes the folloiwing subset of `sourmash gather` metrics:

| Metric | Description |
|--------|-------------|
| Match rank | Order of matches ranked by quality (0 = best match) |
| Species call | Taxonomic classification of the reference genome |
| Query containment ANI | Average Nucleotide Identity for the query intersect section |
| Intersect (bp) | Number of base pairs that overlap between query and reference |
| Remaining (bp) | Base pairs in query not explained by previous matches |
| Fraction original query | Proportion of the original query sequence that matches this reference |
| Fraction unique | Proportion of query that uniquely matches this reference (not shared with others) |
| Reference Genome | Identifier and name of the matching reference genome in database |

The full set of `sourmash gather` metrics for all matches with intersect above 25000bp threshold is output as a `{alias}_sourmash_taxonomy.csv` file.

The workflow uses [GTDB R226](https://sourmash.readthedocs.io/en/latest/databases-md/gtdb226.html) bacterial and [NCBI Funal reference](https://sourmash.readthedocs.io/en/latest/databases-md/ncbi_euks_2025_01.html) genome databases.
Sourmash behavior can be further modified by providing an exclusion list via `--sourmash_db_exclude_list` or by modifying the existing exclusion list present in `data\sourmash_db_exclude_list.txt`.
This parameter accepts a text file containing assembly identifiers (GCA_XXXXXXXXX.V or GCF_XXXXXXXXX.V format) to exclude from the reference databases. 
Excluded genomes are filtered out before species identification, allowing for removal of problematic reference sequences.

#### iii. Antimicrobial resistance (AMR) calling

[ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) is used to detect the presence of acquired resistance genes in the genomes of all bacterial species.
In addition, a subset of well-characterised species/genera will be further analysed using integrated PointFinder to identify specific resistance-associated SNVs.
The following species/genera are included in PointFinder analysis:

* _Campylobacter spp._ 
* _Enterococcus faecalis_
* _Enterococcus faecium_
* _Escherichia coli_
* _Helicobacter pylori_
* _Klebsiella spp._
* _Mycobacterium tuberculosis_
* _Neisseria gonorrhoeae_
* _Salmonella spp._
* _Staphylococcus aureus_
* _Aspergillus fumigatus_
* _Candida albicans_
* _Candidozyma auris_
* _Nakaseomyces glabratus_
* _Saccharomyces cerevisiae_

For purposes of PointFinder analysis, species/genera for a given assembly will be determined by the results of the MLST or Sourmash analysis.
PointFinder analysis will be selected automatically if applicable for a given species call.
PointFinder behaviour can be additionally modified using `--pointfinder_ignore_indels` and `--pointfinder_ignore_stop_codons` parameters.
The role of these options is further described in the [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) documentation.

Optionally, users can replace the default databases with custom ResFinder and PointFinder databases.
Paths to database files can be provided using `--resfinder_db` and `--pointfinder_db`.
Both databases should be KMA indexed as described in the [ResFinder documentation](https://bitbucket.org/genomicepidemiology/resfinder/src/master/).
Additionally, when using a database with added PointFinder species, users need to provide a CSV file mapping Sourmash species names to PointFinder species names.
The file path can be specified using `--sourmash_pointfinder_mapping`.
This CSV file must contain two columns: 'Sourmash species' and 'PointFinder species'.

For example:

```
Sourmash species,PointFinder species
Klebsiella pneumoniae,klebsiella
Klebsiella aerogenes,klebsiella
Candida albicans,candida_albicans
```

In cases where custom databases are not specified, the workflow will use [ResFinder database](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) v2.6.0, [DisinFinder database](https://bitbucket.org/genomicepidemiology/disinfinder_db/src/master/) v2.0.1, and a custom [PointFinder database](https://bitbucket.org/genomicepidemiology/pointfinder_db/src/master/) that has been extended to include fungal mutations.
Fungal AMR detection has been implemented for five key species by using single nucleotide SNVs from [FungAMR DB](https://card.mcmaster.ca/fungamrhome).


#### iv. _Salmonella_ serotyping
Samples identified as salmonella from the MLST step will undergo serotyping and antigenic profile prediction analysis using [SeqSero2](https://github.com/denglab/SeqSero2). The analysis will give predictions on:

* Serotype
* Antigenic profile
* Sub-species identification
* O antigen
* H1 antigen (_fliC_)
* H2 antigen (_fljB_)
