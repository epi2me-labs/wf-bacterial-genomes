<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->
### 1. Concatenates input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities.

### 2a. De-novo assembly

#### i. Assembly

[Flye](https://github.com/fenderglass/Flye) is used to create a draft assembly from the FASTQ reads. This will run by default on the `--nano-raw` paramter for flye. Additional configuration can be performed using `--flye_opts` parameter. 

#### ii. Polishing

The draft assembly from flye is then polished using [medaka](https://github.com/nanoporetech/medaka). This step will attempt to correct any errors that were introduced during the de-novo assembly process. 

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

### 3. Annotations

Regions of interest within your assembly are identified and annotated using [Prokka](https://github.com/tseemann/prokka). By default, prokka will run with it's default databases, but users can refine the annotation using the `--prokka_opts` command. **NOTE** The workflow does not current accept any additional files sent to prokka such as GBK or GFF files.

### 4. Isolates mode (optional)

#### i. Multi-locus sequence typing (MLST)

MLST is a common technique used to help characterise your bacterial isolate, by using allelic variation from internal DNA fragments of 6-7 house keeping genes. Typing schemes for specific species and genera are found on [PubMLST](https://pubmlst.org/) and are pre-loaded into this workflow. [MLST](https://github.com/tseemann/mlst) will try to infer the correct typing scheme to use by scanning the assembly and subsequently identify the allele variant found.

#### ii. Antimicrobial resistance (AMR) calling

[ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) is used to identify genes/SNVs associated with AMR in your assembly. Assemblies of any species will be searched for the detection of acquired resistance genes, however SNVs conferring resistance are only available to a few well characterised species/genera. These are:
* Campylobacter spp.
* Enterococcus faecalis
* Enterococcus faecium
* Escherichia coli
* Helicobacter pylori
* Klebsiella spp. 
* Mycobacterium tuberculosis
* Neisseria gonorrhoeae
* Salmonella spp.
* Staphylococcus aureus

The species/genera of your assembly will be detected from the results of the MLST step and SNV will be selected automatically if applicable.


#### iii. Salmonella serotyping
Samples identified as salmonella from the MLST step will undergo serotyping and antigenic profile prediction analysis using [SeqSero2](https://github.com/denglab/SeqSero2). The analysis will give predictions on:
* Serotype
* Antigenic profile
* Sub-species identification
* O antigen
* H1 antigen (fliC)
* H2 antigen (fljB)

