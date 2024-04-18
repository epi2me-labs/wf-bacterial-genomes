"""Process resfinder results."""

import os
import re

from Bio.Blast import NCBIXML
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def process_resfinder(resfinder_data, amr_type="resfinder"):
    """Select right columns and split loc column."""
    if len(resfinder_data) == 0:
        resfinder_data = pd.DataFrame(columns=[
            'Contig', 'Resistance gene', 'Start', 'End', 'Phenotype',
            'Identity/Nucleotide',
            'Coverage/AA', 'Accession no./PMID', 'Source'])
        return resfinder_data
    loc_split = resfinder_data['Position in contig'].str.split(
        "\\.\\.", n=1, expand=True)
    resfinder_data = resfinder_data.assign(
        Start=loc_split[0], End=loc_split[1])
    resfinder_data = resfinder_data.loc[:, (
        'Contig', 'Resistance gene', 'Start', 'End', 'Phenotype', 'Identity',
        'Coverage', 'Accession no.')]
    resfinder_data = resfinder_data.assign(Type=amr_type)
    out_columns = [
        'Contig', 'Resistance gene', 'Start', 'End', 'Phenotype',
        'Identity/Nucleotide',
        'Coverage/AA', 'Accession no./PMID', 'Source']
    resfinder_data = resfinder_data.set_axis(out_columns, axis=1, copy=False)
    return resfinder_data


# xml_file="pointfinder_blast/tmp/out_pncA-promoter-size-107bp.xml"
# Get the position of a mutation in a sequence
def read_pointfinder_xml(xml_file):
    """Get the position of a mutation in the input contig."""
    records_p = NCBIXML.parse(open(xml_file, "r"))
    results = {
        "query_start": [], "query_end": [], "subject_start": [],
        "subject_end": [], "query": []}

    for record in records_p:
        if record.alignments:
            for align in record.alignments:
                for hsp in align.hsps:
                    results['query_start'].append(hsp.query_start)
                    results['query_end'].append(hsp.query_end)
                    results['subject_start'].append(hsp.sbjct_start)
                    results['subject_end'].append(hsp.sbjct_end)
                    results['query'].append(record.query)

    results_df = pd.DataFrame(results)
    # Switch start and end if start is greater than end
    if (results_df.loc[0, 'query_start'] > results_df.loc[0, 'query_end']):
        results_df['query_start'], results_df['query_end'] = results_df[
            'query_end'], results_df['query_start']
        results_df['subject_start'], results_df['subject_end'] = results_df[
            'subject_end'], results_df['subject_start']
    return results_df


def extract_point_bp(subject, mutation):
    """Extract mutation from pointfinder output."""
    # subject = "pncA-promoter-size-107bp.xml"
    # resfinder 4.4.3 size not in table so subject is now xml
    # mutation = "p.H57D"
    mutation_bp = int(re.search(r'\d+', mutation).group())*3
    if ("promoter" in subject):
        prom_size = int(re.search(r"(\d+)bp.xml$", subject).groups()[0])
        # The noted mutation is from the start of the CDS, if a promoter is
        # present then this should take that into account
        final_bp = prom_size + mutation_bp
        return (final_bp)
    else:
        return (mutation_bp)


def get_global_mutation_position(q_start, q_end, s_start, s_end, loc_pos):
    """Get global mutation position on query depending on query strand orientation."""
    if s_start < s_end:
        glob_loc_start = q_start + loc_pos
        glob_loc_end = q_start + loc_pos + 2
    else:
        glob_loc_start = q_end - loc_pos
        glob_loc_end = q_end - loc_pos - 2
    return (glob_loc_start, glob_loc_end)


# Maybe add a test for this?
def convert_pointfinder_row(database_location, row):
    """Convert point finder row."""
    # construct path to xml file
    (xml_fname,) = [
        f
        for f in os.listdir(database_location)
        if f.startswith(f"out_{row['Sequence']}")
    ]
    xml_path = os.path.join(database_location, xml_fname)
    pointfinder_blast_data = read_pointfinder_xml(xml_path)
    pointfinder_blast_data['Subject'] = row['Sequence']
    pointfinder_blast_data['Mutation'] = row['Mutation']
    # TODO: This works for mutations in CDS, but not sure how it will work
    # with mutations in promoter regions

    pointfinder_blast_data['Local_loc_bp'] = pointfinder_blast_data.apply(
        lambda row: extract_point_bp(xml_path, row["Mutation"]), axis=1
        )
    # Find start and end of mutation on query depending on orientation of gene
    pointfinder_blast_data[
        ['Global_loc_bp', 'Global_loc_bp_end']
        ] = pointfinder_blast_data.apply(
            lambda row: get_global_mutation_position(
                row["query_start"],
                row["query_end"],
                row["subject_start"],
                row["subject_end"],
                row["Local_loc_bp"]
            ), axis=1, result_type="expand"
        )
    return pointfinder_blast_data


def process_pointfinder(pointfinder_data, database_location):
    """Select right columns and split loc column."""
    # TODO: This is all messy and needs refactoring
    if len(pointfinder_data) == 0:
        pointfinder_data = pd.DataFrame(columns=[
            'Contig', 'Start', 'End', 'Phenotype', 'Identity/Nucleotide',
            'Coverage/AA', 'Accession no./PMID', 'Source'])
        return pointfinder_data
    try:
        gene_name = pointfinder_data['Mutation'].str.split('-').str[0]
    except IndexError:
        gene_name = pointfinder_data['Mutation']
    pointfinder_data["Resistance gene"] = gene_name
    loc_split = pointfinder_data['Mutation'].str.split(
        " ", n=1, expand=True)
    pointfinder_data = pointfinder_data.assign(
        Sequence=loc_split[0], Mutation=loc_split[1])
    pointfinder_data['Query_loc'] = 0
    df_list = []
    # multiple entries can be found in xml so iterate over rows
    for _, row in pointfinder_data.iterrows():
        pointfinder_blast_data = convert_pointfinder_row(
            database_location, row)
        df_list.append(pointfinder_blast_data)

    concat_df = pd.concat(df_list)
    concat_df_merge = concat_df.merge(pointfinder_data, on='Mutation')
    concat_df_merge_tidy = concat_df_merge.iloc[
        :, [4, 8, 9, 12, 10, 11, 13, 14]]
    concat_df_merge_tidy = concat_df_merge_tidy.assign(Type='pointfinder')
    out_columns = [
        'Contig', 'Start', 'End', 'Phenotype', 'Identity/Nucleotide',
        'Coverage/AA', 'Accession no./PMID', 'Resistance gene', 'Source']
    concat_df_merge_tidy = concat_df_merge_tidy.set_axis(
        out_columns, axis=1, copy=False)
    return concat_df_merge_tidy


def main(args):
    """Run entry point."""
    logger = get_named_logger("process_amr")
    resfinder_data = pd.read_csv(args.resfinder_file, sep="\t")
    resfinder_results = process_resfinder(resfinder_data)
    if (args.pointfinder_file != "NOT_RUN"):
        pointfinder_data = pd.read_csv(args.pointfinder_file, sep="\t")
        # Drop duplicates - entries will be recovered in process_pointfinder
        pointfinder_data = pointfinder_data.drop_duplicates()
        pointfinder_results = process_pointfinder(
            pointfinder_data, args.database_location)
    else:
        pointfinder_results = pd.DataFrame()
    if (args.disinf_file != "NOT_RUN"):
        disinf_data = pd.read_csv(args.disinf_file, sep="\t")
        disinf_results = process_resfinder(
            disinf_data, amr_type="Disinfectant")
    else:
        disinf_results = pd.DataFrame()
    amr_results = pd.concat([
        resfinder_results, pointfinder_results, disinf_results])
    amr_results.to_csv(args.output, sep="\t", index=False, na_rep='N/A')
    logger.info(f"Written amr-results to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("process_amr")
    parser.add_argument(
        "--resfinder_file",
        help="Resfinder tab delimited results file.")
    parser.add_argument(
        "--disinf_file", default="NOT_RUN",
        help="Disinfectant resistance finder tab delimited results file.")
    parser.add_argument(
        "--pointfinder_file", default="NOT_RUN",
        help="Pointfinder tab delimited results file.")
    parser.add_argument(
        "--database_location",
        help="Location of pointfinder results tmp storage.")
    parser.add_argument(
        "--species_id",
        help="Species id for pointfinder database"
        "(e.g mycobacterium_tuberculosis).")
    parser.add_argument(
        "--output", default="Results_tab_processed.txt",
        help="Destination file for processed resfinder results.")

    return parser
