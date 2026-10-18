# Running wf-bacterial-genomes with funannotate (fungal genomes)

This fork adds [funannotate](https://github.com/nextgenusfs/funannotate) as an
alternative to Bakta for the **Annotation** step. IT can be used for **fungal**
assemblies; Bakta remains the default for bacteria.

```
reads -> Flye -> dnaapler -> mini_align -> medaka  -> consensus assembly
                                                          |
                              +---------------------------+-----------------+
                              | --run_bakta (default)     | --run_funannotate
                              v                           v
                            Bakta              funannotate mask -> predict -> annotate
                              \                          /
                               +--> GFF3 -> report "Annotations" table
```

Notes: piu piu funannotate is run **without RNA-seq evidence** (ab initio + protein evidence
from the database, Augustus/SNAP/GlimmerHMM; GeneMark is not available because
of its licence for edu only).

## 1. Requirements

| Requirement | Notes |
| --- | --- |
| Nextflow >= 23.04.2 and Java | tested with 24.10.5 / 25.10.5 |
| Docker | the funannotate image is built locally (see 2) |
| Host resources | >= 8 CPUs / 32 GB RAM (16 CPUs / 64 GB recommended); the funannotate process requests `--threads` CPUs and 16 GB |
| funannotate database | built beforehand manually, because go.obo did not work with `funannotate setup` (see 3) |

## 2. Build the container image (once)

The image is built from the pinned conda environment
[`funannotate-env/funannotate.yml`](funannotate-env/funannotate.yml) (do not
edit or regenerate it).

```bash
cd wf-bacterial-genomes
docker build --platform linux/amd64 \
    -t funannotate-wf:1.8.17-2 \
    -f funannotate-env/Dockerfile funannotate-env/
```

The tag must match `params.wf.container_tag_funannotate` in
[`nextflow.config`](nextflow.config). **Never reuse a tag**: if you change the
Dockerfile or the yml, build `1.8.17-3` and update `nextflow.config`.

Sanity check (should only miss GeneMark, `emapper.py` and `signalp`, which are
optional/licensed):

```bash
docker run --rm --user $(id -u):$(id -g) funannotate-wf:1.8.17-2 \
    funannotate check --show-versions
```

> The image is local to the machine that built it. It will **not** work with
> Singularity, other hosts or AWS Batch until it is pushed to a registry and the
> `withLabel:funannotate` container in `nextflow.config` is changed to the
> registry name.

## 3. Database

```bash
export FUNANNOTATE_DB=/data/datadrive/5.IGF/Resources/Databases/funannotate_db
```

The database is **not** baked into the image; it is passed with
`--funannotate_db` and exported as `FUNANNOTATE_DB` inside the task.

* It contains BUSCO lineage `dikarya` (the default used by `funannotate predict`).
  For any other lineage add e.g. `--funannotate_opts "--busco_db fungi"` and make
  sure that lineage exists in the database.
* `funannotate setup` writes **absolute paths** into the database
  (`funannotate-db-info.txt`, `trained_species/*/info.json`). Nextflow mounts the
  database at its real host path, so this works as-is, but if you move or copy
  the database you must rebuild/fix it, and any manual `docker run` must mount it
  at the **same path** (`-v $DB:$DB:ro`), not at e.g. `/db`.

## 4. Run

Fungal reads are given like any other run of this workflow (a directory of
FASTQ, or a sample sheet with `barcode,alias,type`).

```bash
nextflow run . -c "$PWD/nextflow.config" \
    --fastq /path/to/fastq_pass \
    --sample_sheet /path/to/sample_sheet.csv \
    --run_bakta false \
    --run_funannotate \
    --funannotate_species "Aspergillus fumigatus" \
    --funannotate_db /data/datadrive/5.IGF/Resources/Databases/funannotate_db \
    --isolates false \
    --run_plasmid_id false \
    --threads 16 \
    -profile standard
```

Preview what would be scheduled without running anything:

```bash
nextflow run . -preview --fastq test_data/fungi \
    --run_bakta false --run_funannotate \
    --funannotate_species "Aspergillus fumigatus" \
    --funannotate_db $FUNANNOTATE_DB -profile standard
```

Notes:

* `--run_bakta false` is required: if both `--run_bakta` and `--run_funannotate`
  are true, **Bakta wins**.
* `--isolates false --run_plasmid_id false` switch off the bacteria-specific
  characterisation (MLST, Sourmash species ID, ResFinder, SeqSero2, MOB-suite),
  which is meaningless for fungi. They are on by default in the workflow.
* Use a species close to your organism for `--funannotate_species`; if it matches
  a pre-trained Augustus species in `trained_species/` (e.g. `aspergillus_fumigatus`)
  that training is used, otherwise Augustus is trained from BUSCO.
* Contig names must be short (funannotate default limit is 16 characters). Flye's
  `contig_N` names are fine.

### Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--run_funannotate` | `false` | Run funannotate `mask` + `predict` + `annotate` on the consensus assembly. |
| `--funannotate_species` | `null` | **Required.** Species name given to `funannotate -s`, e.g. `"Aspergillus fumigatus"`. |
| `--funannotate_db` | `null` | **Required.** Path to the funannotate database. |
| `--funannotate_opts` | `null` | Extra arguments for `funannotate predict`, e.g. `"--busco_db fungi"`. |
| `--run_bakta` | `true` | Must be `false` to use funannotate. |
| `--threads` | `3` | CPUs for funannotate (and other processes). |

The workflow stops immediately with `--funannotate_species is required ...` or
`--funannotate_db is required ...` if either is missing.

## 5. Outputs

In `output/` (per sample alias):

| File | Description |
| --- | --- |
| `<alias>.funannotate.gff3` | Final annotated gene models (`funannotate annotate`) |
| `<alias>.funannotate.gbk` | Same annotation in GenBank format |
| `wf-bacterial-genomes-report.html` | Report; the **Annotations** section says "annotated with funannotate" |

In the report table each row is one `mRNA`, `tRNA`, `rRNA` or `ncRNA` feature
(gene/exon/CDS rows are not listed). The gene symbol comes from the parent `gene`
row when present, and the EC number from `EC_number` or `db_xref=EC:...`.

If funannotate fails for a sample, the run continues, a
`funannotate failed for sample '<alias>': <last 20 stderr lines>` warning is
printed, and the report shows *No annotation data available for this sample.*

## 6. Troubleshooting

| Symptom | Cause / fix |
| --- | --- |
| `Found unexpected parameters` | Parameter not declared in `nextflow_schema.json` |
| `Can't find Repeat Database at .../repeats.dmnd` | `--funannotate_db` is wrong or the database was not built with `funannotate setup` |
| `FileNotFoundError ... trained_species/<sp>/augustus` | The database was moved after `funannotate setup`; its stored absolute paths are stale (see 3) |
| `busco database is not found` | The BUSCO lineage is missing from the database, or `--busco_db` names another lineage |
| Annotation table empty but GFF3 exists | Run was started with `-resume` after editing Python in `bin/`; rerun without `-resume` |
| Bakta image still pulled | Only happens when `--run_funannotate` is not set or `--run_bakta` is true |
| Nextflow ignores the repo `nextflow.config` / `out_dir undefined` | `NXF_CONFIG_FILE` points to a missing file; `unset NXF_CONFIG_FILE` |
| `no such image: funannotate-wf:<tag>` | Image not built on this host, or tag differs from `container_tag_funannotate` |
