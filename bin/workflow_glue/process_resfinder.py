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
                    results['query'].append(record.query[:100])

    results_df = pd.DataFrame(results)
    # Switch start and end if start is greater than end
    if (results_df.loc[0, 'query_start'] < results_df.loc[0, 'query_end']):
        results_df['query_start'], results_df['query_end'] = results_df[
            'query_end'], results_df['query_start']
        results_df['subject_start'], results_df['subject_end'] = results_df[
            'subject_end'], results_df['subject_start']
    return results_df


def extract_point_bp(subject, mutation):
    """Extract mutation from pointfinder output."""
    # subject = "pncA-promoter-size-107bp"
    # mutation = "p.H57D"
    mutation_bp = int(re.search(r'\d+', mutation).group())*3
    if ("promoter" in subject):
        bp_re = re.compile("[0-9]*bp")
        prom_split = subject.split("-")
        prom_size = int([x for x in prom_split if bp_re.match(
            x) is not None][0].replace("bp", ""))
        # The noted mutation is from the start of the CDS, if a promoter is
        # present then this should take that into account
        final_bp = prom_size + mutation_bp
        return (final_bp)
    else:
        return (mutation_bp)


# Maybe add a test for this?
def convert_pointfinder_row(database_location, row):
    """Convert point finder row."""
    # construct path to xml file
    xml_path = os.path.join(
        database_location, "out_"+row['Sequence'] + ".xml")
    pointfinder_blast_data = read_pointfinder_xml(xml_path)
    pointfinder_blast_data['Subject'] = row['Sequence']
    pointfinder_blast_data['Mutation'] = row['Mutation']
    # TODO: This works for mutations in CDS, but not sure how it will work
    # with mutations in promoter regions
    bp_position = extract_point_bp(
        pointfinder_blast_data.loc[0, 'Subject'],
        pointfinder_blast_data.loc[0, 'Mutation'])
    pointfinder_blast_data['Local_loc_bp'] = bp_position
    start_loc = pointfinder_blast_data.loc[0, 'query_start']
    end_loc = pointfinder_blast_data.loc[0, 'query_end']
    # Getting round the fact that the res-gene may match in the reverse
    # orientation, and therefore the start and end positions are reversed
    if (start_loc > end_loc):
        pointfinder_blast_data['Global_loc_bp'] = start_loc - bp_position
        # End of codon
        pointfinder_blast_data['Global_loc_bp_end'] = start_loc - \
            bp_position + 2
    else:
        pointfinder_blast_data['Global_loc_bp'] = end_loc + bp_position
        # TODO: should this be a + or -
        pointfinder_blast_data['Global_loc_bp_end'] = start_loc - \
            bp_position + 2
    return pointfinder_blast_data


def process_pointfinder(pointfinder_data, output_file, database_location):
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
    mutation_names = pointfinder_data['Mutation']
    loc_split = pointfinder_data['Mutation'].str.split(
        " ", n=1, expand=True)
    pointfinder_data = pointfinder_data.assign(
        Sequence=loc_split[0], Mutation=loc_split[1])
    pointfinder_data['Query_loc'] = 0
    df_list = []
    for key, row in pointfinder_data.iterrows():
        pointfinder_blast_data = convert_pointfinder_row(
            database_location, row)
        df_list.append(pointfinder_blast_data)

    concat_df = pd.concat(df_list)
    # split the mutation and get the base position
    concat_df_merge = concat_df.merge(pointfinder_data, on='Mutation')
    concat_df_merge_tidy = concat_df_merge.iloc[
        :, [4, 8, 9, 12, 10, 11, 13, 14]]
    concat_df_merge_tidy.index = mutation_names
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
        pointfinder_results = process_pointfinder(
            pointfinder_data, args.output, args.database_location)
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
