"""Collect results into wf.Sample."""

import json
import os

import pandas as pd
import vcf
from workflow_glue.models.common import SampleType
from workflow_glue.models.custom import (
    Annotation, AntimicrobialResistance, Assembly,
    Contig, ContigType, Coverage, FastqStats, FlyeContigStats, MLST,
    PlasmidID, PlasmidIdentification, ResultsContents,
    Sample, SequenceTypeSchema, Serotype, SourmashExcludedAssemblies,
    SourmashMatch, SpeciesIdentification, Variant
)

from .parsers import (  # noqa: ABS101
    parse_bakta_gff3,
    parse_generic_gff3
)
from .process_resfinder_iso import (  # noqa: ABS101
    parse_resfinder
)
from .util import get_named_logger, wf_parser  # noqa: ABS101


def gather_sample_files(alias, data_dir):
    """Collect files required for the model."""
    files = {
        "depth": os.path.join(
            data_dir, f"{alias}.total.regions.bed.gz"),
        "variants": os.path.join(
            data_dir, f"{alias}.medaka.vcf.gz"),
        "bakta": os.path.join(
            data_dir, f"{alias}.bakta.gff3"),
        "funannotate": os.path.join(
            data_dir, f"{alias}.funannotate.gff3"),
        "mlst": os.path.join(
            data_dir, f"{alias}.mlst.json"),
        "amr": os.path.join(
            data_dir, f"{alias}_resfinder_results/{alias}_resfinder.json"),
        "pointfinder_genes": os.path.join(
            data_dir, f"{alias}_resfinder_results/pointfinder_genes_analysed.txt"),
        "disinfinder_table": os.path.join(
            data_dir, f"{alias}_resfinder_results/DisinFinder_results_table.txt"),
        "fastcat": os.path.join(
            data_dir, "fastcat_stats/per-read-stats.tsv.gz"),
        "flye_stats": os.path.join(
            data_dir, f"{alias}.flye_stats.tsv"),
        "serotype": os.path.join(
            data_dir, f"{alias}.serotype_results.tsv"),
        "taxonomy": os.path.join(
            data_dir, f"{alias}_sourmash_taxonomy.csv"),
        "excluded_assemblies": os.path.join(
            data_dir, "sourmash_picklist_excluded.txt"),
        "mobsuite_contig_report": os.path.join(
            data_dir, f"{alias}_mob_results/contig_report.txt"),
        "mobsuite_mobtyper": os.path.join(
            data_dir, f"{alias}_mob_results/mobtyper_results.txt")
    }

    # Return none if file does not exist
    files = {
        section: (
            file if os.path.exists(file) else None
            ) for section, file in files.items()
    }

    return files


def fastcat_stats(per_read_file):
    """Parse per-read-statistics.tsv fastcat output."""
    if not per_read_file or not os.path.exists(per_read_file):
        return FastqStats()
    df = pd.read_csv(per_read_file, sep="\t", compression='gzip')
    if df.empty:
        return FastqStats()
    n_seqs = len(df)
    n_bases = int(df['read_length'].sum())
    min_length = int(df['read_length'].min())
    max_length = int(df['read_length'].max())
    median_length = float(df['read_length'].median())
    mean_quality = float(df['mean_quality'].mean())
    median_quality = float(df['mean_quality'].median())

    return FastqStats(
        n_seqs=n_seqs,
        n_bases=n_bases,
        min_length=min_length,
        max_length=max_length,
        median_length=median_length,
        mean_quality=mean_quality,
        median_quality=median_quality
    )


def extract_flye_field(row, key, cast_type=str):
    """Get fields from flye data."""
    if key in row and pd.notna(row[key]):
        return cast_type(row[key])
    return None


def parse_flye_stats(flye_stats_file):
    """Parse Flye assembly statistics TSV file."""
    df = pd.read_csv(flye_stats_file, sep="\t")
    flye_contigs = []
    for _, row in df.iterrows():
        flye_contigs.append(
            FlyeContigStats(
                seq_name=extract_flye_field(row, "#seq_name", str),
                length=extract_flye_field(row, "length", int),
                coverage=extract_flye_field(row, "cov.", float),
                circular=extract_flye_field(row, "circ.", str),
                repeat=extract_flye_field(row, "repeat", str),
                multiplicity=extract_flye_field(row, "mult.", int),
            )
        )
    return flye_contigs


def parse_mlst(mlst_file):
    """Extract schema and sequence type from MLST json."""
    with open(mlst_file, "r") as f:
        mlst = json.loads(f.read())[0]

    alleles = []
    if mlst["alleles"]:
        for schema_id, allele in mlst["alleles"].items():
            alleles.append(SequenceTypeSchema(
                schema_identifier=schema_id,
                allele_variant=allele
            ))
    result = MLST(
        detected_species=mlst["scheme"],
        sequence_type=mlst["sequence_type"],
        typing_schema=alleles
    )

    return result


def parse_sourmash_excluded_genomes(excluded_file):
    """Parse Sourmash excluded assemblies file."""
    if not excluded_file or not os.path.exists(excluded_file):
        return SourmashExcludedAssemblies(excluded_assemblies=[], total_excluded=0)

    excluded_assemblies = []
    with open(excluded_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                excluded_assemblies.append(line)

    return SourmashExcludedAssemblies(
        excluded_assemblies=excluded_assemblies,
        total_excluded=len(excluded_assemblies)
    )


def extract_species_from_lineage(lineage_str):
    """Extract species name from taxonomic lineage string."""
    if lineage_str is None or pd.isna(lineage_str):
        return "Unknown"
    if isinstance(lineage_str, str) and not lineage_str.strip():
        return "Unknown"
    # Extract species
    fields = str(lineage_str).split(';')
    genus = None
    for field in reversed(fields):
        f = field.strip()
        if f.startswith('s__'):  # avoids matching 'ss__'
            species = f[3:].strip()
            if species:
                return species  # checks for empty species string
        elif f.startswith('g__'):
            genus = f[3:].strip()
    return genus or "Unknown"


def parse_sourmash_taxonomy(taxonomy_file):
    """Parse sourmash annotated results."""
    df = pd.read_csv(taxonomy_file)
    if df.empty:
        return SpeciesIdentification(detected_matches=[])
    matches = []
    # Maximum of 5 top hits are returned
    for _, row in df.head(5).iterrows():
        lineage = str(row.get('lineage', '') or '')
        match = SourmashMatch(
            rank=int(row.get('gather_result_rank', len(matches) + 1)),
            reference_name=str(row.get('name', 'Unknown')),
            species_call=extract_species_from_lineage(lineage),
            intersect_bp=(
                int(row['intersect_bp'])
                if pd.notna(row.get('intersect_bp')) else None
            ),
            f_orig_query=(
                float(row['f_orig_query'])
                if pd.notna(row.get('f_orig_query')) else None
            ),
            f_unique_to_query=(
                float(row['f_unique_to_query'])
                if pd.notna(row.get('f_unique_to_query')) else None
            ),
            query_containment_ani=(
                float(row['query_containment_ani'])
                if pd.notna(row.get('query_containment_ani')) else None
            ),
            remaining_bp=(
                int(row['remaining_bp'])
                if pd.notna(row.get('remaining_bp')) else None
            )
        )
        matches.append(match)
    return SpeciesIdentification(detected_matches=matches)


def parse_serotyping(serotype_file):
    """Extract serotyping information from SeqSero2."""
    sero_df = pd.read_csv(serotype_file, sep="\t")
    # columns always present in seqsero output
    serotype = Serotype(
        predicted_serotype=str(
            sero_df["Predicted serotype"].squeeze()),
        predicted_antigenic_profile=str(
            sero_df["Predicted antigenic profile"].squeeze()),
        predicted_identification=str(
            sero_df["Predicted identification"].squeeze()),
        o_antigen_prediction=str(
            sero_df["O antigen prediction"].squeeze()),
        h1_antigen_prediction=str(
            sero_df["H1 antigen prediction(fliC)"].squeeze()),
        h2_antigen_prediction=str(
            sero_df["H2 antigen prediction(fljB)"].squeeze()),
        note=str(
            sero_df["Note"].squeeze())
    )
    return serotype


def annotation_stats(gff_file, tool="bakta"):
    """Extract annotation information from a Bakta or funannotate GFF3 file."""
    if tool == "bakta":
        annot_df = parse_bakta_gff3(gff_file)
    else:
        annot_df = parse_generic_gff3(
            gff_file, keep_types={"mRNA", "tRNA", "rRNA", "ncRNA"})
    annot_df = annot_df.rename(columns={"EC number": "ec_number"})

    annotations = []
    for _, row in annot_df.iterrows():
        annotations.append(Annotation(
            contig=str(row["contig"]) if str(row["contig"]) != "-" else None,
            ID=str(row["ID"]) if str(row["ID"]) != "-" else None,
            start=int(str(row["start"])) if str(row["start"]) != "-" else None,
            end=int(str(row["end"])) if str(row["end"]) != "-" else None,
            strand=str(row["strand"]) if str(row["strand"]) != "-" else None,
            gene=str(row["gene"]) if str(row["gene"]) != "-" else None,
            product=str(row["product"]) if str(row["product"]) != "-" else None,
            ec_number=str(row["ec_number"]) if str(row["ec_number"]) != "-" else None
        ))
    return annotations


def contig_stats(total_coverage):
    """Get coverage details for assembly."""
    contigs = []
    depth_df = pd.read_csv(
        total_coverage, sep="\t", names=["ref", "start", "end", "depth"]
    )
    for contig, df in depth_df.groupby("ref"):
        coverage = Coverage(
            counts=int(df["depth"].sum()),
            median=float(df["depth"].median()),
            mean=float(df["depth"].mean()),
            minimum=int(df["depth"].min()),
            maximum=int(df["depth"].max())
        )
        contigs.append(Contig(
            name=str(contig),
            length=int(df["end"].max()),
            coverage=coverage
        ))
    return contigs


def variant_stats(vcf_file):
    """Extract basic variant information from vcf file."""
    vcf_reader = vcf.Reader(filename=vcf_file)
    variants = []
    for record in vcf_reader:
        for alt in record.ALT:
            variants.append(Variant(
                contig=record.CHROM,
                pos=record.POS,
                ref=record.REF,
                alt=str(alt),
                depth=record.INFO["DP"]
            ))
    return variants


def assembly_stats(params_data, files):
    """Gather assembly statistics for sample."""
    contigs = []
    if files["depth"]:
        contigs = contig_stats(files["depth"])

    variants = []
    if files["variants"]:
        variants = variant_stats(files["variants"])

    annotations = []
    if params_data.get("run_bakta", False) and files["bakta"]:
        annotations = annotation_stats(files["bakta"])
    elif params_data.get("run_funannotate", False) and files["funannotate"]:
        annotations = annotation_stats(files["funannotate"], tool="funannotate")

    flye_stats = None
    if files["flye_stats"]:
        flye_stats = parse_flye_stats(files["flye_stats"])

    assembly = Assembly(
        reference=params_data.get("reference"),
        annotations=annotations,
        contig=contigs,
        variants=variants,
        flye_stats=flye_stats
    )

    return assembly


def antimicrobial_stats(
    resfinder_json, pointfinder_genes_file=None, disinfinder_table_file=None
):
    """Gather amr results from resfinder."""
    with open(resfinder_json) as f:
        resfinder_data = json.loads(f.read())

    parsed_data = parse_resfinder(resfinder_data, disinfinder_table_file)

    if pointfinder_genes_file:
        pointfinder_available_genes = parse_pointfinder_genes(pointfinder_genes_file)
    else:
        pointfinder_available_genes = None

    results = AntimicrobialResistance(
        databases=parsed_data["databases"],
        seq_regions=parsed_data["seq_regions"],
        seq_variations=parsed_data["seq_variations"],
        phenotypes=parsed_data["phenotypes"],
        provided_species=parsed_data["provided_species"],
        pointfinder_available_genes=pointfinder_available_genes
    )
    return results


def extract_mobsuite_field(row, key, cast_type):
    """Get fields from MOB-suite data, handling NaN and dash placeholders."""
    if key in row and pd.notna(row[key]):
        value = str(row[key]).strip()
        if value != '-' and value != '':
            return cast_type(value)
    return None


def parse_mobtyper_results(mobtyper_file):
    """Parse MOB-suite mobtyper results."""
    df = pd.read_csv(mobtyper_file, sep="\t")

    mobtyper_data = {}
    for _, row in df.iterrows():
        sample_id = extract_mobsuite_field(row, 'sample_id', str)
        if sample_id and ':' in sample_id:
            cluster_id = sample_id.split(':')[1]
            mobtyper_data[cluster_id] = {
                'mpf_type':
                    extract_mobsuite_field(row, 'mpf_type', str),
                'predicted_mobility':
                    extract_mobsuite_field(row, 'predicted_mobility', str),
                'predicted_host_range_rank':
                    extract_mobsuite_field(
                        row, 'predicted_host_range_overall_rank', str),
                'predicted_host_range_name':
                    extract_mobsuite_field(
                        row, 'predicted_host_range_overall_name', str),
                'observed_host_range_rank':
                    extract_mobsuite_field(row, 'observed_host_range_ncbi_rank', str),
                'observed_host_range_name':
                    extract_mobsuite_field(row, 'observed_host_range_ncbi_name', str)
            }

    return mobtyper_data


def parse_mobsuite_results(mobsuite_file, mobtyper_file=None):
    """Parse MOB-suite contig report."""
    df = pd.read_csv(mobsuite_file, sep="\t")

    if df.empty:
        return PlasmidIdentification()

    mobtyper_data = {}
    if mobtyper_file:
        mobtyper_data = parse_mobtyper_results(mobtyper_file)

    detected_contigs = []
    for _, row in df.iterrows():
        contig_id_raw = extract_mobsuite_field(row, 'contig_id', str)
        if contig_id_raw:
            if '_basecall_model=' in contig_id_raw:
                contig_id = contig_id_raw.split('_basecall_model=')[0]
            else:
                contig_id = contig_id_raw.split(' ')[0]
        else:
            contig_id = None

        molecule_type_str = extract_mobsuite_field(row, 'molecule_type', str)
        molecule_type = None
        if (
            molecule_type_str and
            molecule_type_str.lower() in [e.value for e in ContigType]
        ):
            molecule_type = ContigType(molecule_type_str.lower())

        # Get initial values and primary cluster ID
        mpf_type = extract_mobsuite_field(row, 'mpf_type', str)
        predicted_mobility = extract_mobsuite_field(row, 'predicted_mobility', str)
        primary_cluster_id = extract_mobsuite_field(row, 'primary_cluster_id', str)
        predicted_host_range_rank = None
        predicted_host_range_name = None
        observed_host_range_rank = None
        observed_host_range_name = None

        # mobtyper_results.txt (if present) contains more detailed info for these fields
        if primary_cluster_id in mobtyper_data:
            mobtyper_info = mobtyper_data[primary_cluster_id]
            if mobtyper_info['mpf_type']:
                mpf_type = mobtyper_info['mpf_type']
            if mobtyper_info['predicted_mobility']:
                predicted_mobility = mobtyper_info['predicted_mobility']
            if mobtyper_info.get('predicted_host_range_rank'):
                predicted_host_range_rank = mobtyper_info['predicted_host_range_rank']
            if mobtyper_info.get('predicted_host_range_name'):
                predicted_host_range_name = mobtyper_info['predicted_host_range_name']
            if mobtyper_info.get('observed_host_range_rank'):
                observed_host_range_rank = mobtyper_info['observed_host_range_rank']
            if mobtyper_info.get('observed_host_range_name'):
                observed_host_range_name = mobtyper_info['observed_host_range_name']

        plasmid = PlasmidID(
            contig_id=contig_id,
            contig_type=molecule_type,
            primary_cluster_id=primary_cluster_id,
            size=extract_mobsuite_field(
                row, 'size', int),
            rep_type=extract_mobsuite_field(
                row, 'rep_type(s)', str),
            relaxase_type=extract_mobsuite_field(
                row, 'relaxase_type(s)', str),
            mpf_type=mpf_type,
            orit_type=extract_mobsuite_field(
                row, 'orit_type(s)', str),
            predicted_mobility=predicted_mobility,
            mash_nearest_neighbor=extract_mobsuite_field(
                row, 'mash_nearest_neighbor', str),
            mash_neighbor_distance=extract_mobsuite_field(
                row, 'mash_neighbor_distance', float),
            mash_neighbor_identification=extract_mobsuite_field(
                row, 'mash_neighbor_identification', str),
            predicted_host_range_rank=predicted_host_range_rank,
            predicted_host_range_name=predicted_host_range_name,
            observed_host_range_rank=observed_host_range_rank,
            observed_host_range_name=observed_host_range_name
        )
        detected_contigs.append(plasmid)

    plasmid_count = 0
    chromosome_count = 0
    for contig in detected_contigs:
        if contig.contig_type == ContigType.PLASMID:
            plasmid_count += 1
        elif contig.contig_type == ContigType.CHROMOSOME:
            chromosome_count += 1

    return PlasmidIdentification(
        detected_contigs=detected_contigs,
        n_plasmid_contigs=plasmid_count,
        n_chromosome_contigs=chromosome_count
        )


def parse_pointfinder_genes(pointfinder_genes_file):
    """Parse PointFinder available genes file."""
    if not pointfinder_genes_file or not os.path.exists(pointfinder_genes_file):
        return []
    genes = []
    with open(pointfinder_genes_file, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene:
                genes.append(gene)
    return genes


def main(args):
    """Run the entry point."""
    logger = get_named_logger("collect_results")
    with open(args.params, "r") as f:
        params_data = json.loads(f.read())

    alias = args.alias
    files = gather_sample_files(args.alias, args.data_dir)
    assembly = assembly_stats(params_data, files)

    if files["mlst"]:
        sequence_type = parse_mlst(files["mlst"])
    else:
        sequence_type = MLST()

    if files["taxonomy"]:
        species_identification = parse_sourmash_taxonomy(files["taxonomy"])
    else:
        species_identification = SpeciesIdentification(detected_matches=[])

    if files["amr"]:
        antimicrobial_details = antimicrobial_stats(
            files["amr"],
            files["pointfinder_genes"],
            files["disinfinder_table"]
        )
    else:
        antimicrobial_details = AntimicrobialResistance()

    if files["serotype"]:
        serotype = parse_serotyping(files["serotype"])
    else:
        serotype = Serotype()

    if files["fastcat"]:
        fastcat = fastcat_stats(
            files["fastcat"]
        )
    else:
        fastcat = FastqStats()

    if files["excluded_assemblies"]:
        sourmash_excluded = parse_sourmash_excluded_genomes(
            files["excluded_assemblies"]
        )
    else:
        sourmash_excluded = SourmashExcludedAssemblies(
            excluded_assemblies=[],
            total_excluded=0
        )

    if files["mobsuite_contig_report"]:
        plasmid_identification = parse_mobsuite_results(
            files["mobsuite_contig_report"],
            files["mobsuite_mobtyper"]
        )
    else:
        plasmid_identification = PlasmidIdentification()

    results = ResultsContents(
        antimicrobial_resistance=antimicrobial_details,
        assembly=assembly,
        sequence_typing=sequence_type,
        species_identification=species_identification,
        serotyping=serotype,
        fastq=fastcat,
        sourmash_excluded_genomes=sourmash_excluded,
        plasmid_identification=plasmid_identification
    )

    sample = Sample(
        alias=alias,
        sample_type=SampleType(args.type),
        results=results,
        barcode=args.barcode,
        additional_identifiers=[],
        sample_checks=[],
        config=params_data
    )

    sample.to_json(args.output)
    logger.info(f"results collected and written to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("collect_results")
    parser.add_argument("--output", help="Report output filename")
    parser.add_argument("--params", required=True)
    parser.add_argument("--alias", required=True)
    parser.add_argument("--barcode", required=True)
    parser.add_argument("--data_dir", required=True, help="Analysis results directory")
    parser.add_argument("--type")
    return parser
