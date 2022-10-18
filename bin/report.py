#!/usr/bin/env python
"""Create workflow report."""

import argparse
import base64
import glob
import io
import os
from pathlib import Path


from aplanat import report
from aplanat.components import bcfstats
from aplanat.components import depthcoverage, fastcat
from aplanat.components import simple as scomponents
from aplanat.util import Colors
from Bio import SeqIO
from bokeh.layouts import gridplot, layout
from bokeh.models import Panel, Tabs
from dna_features_viewer import BiopythonTranslator
import pandas as pd


def collate_stats(target_glob: str, input_sep="\t", filename_merge=False,
                  **kwargs) -> pd.DataFrame:
    r"""Collate stats files.

    Args:
        target_glob (str): Target glob to get files to read
        input_sep (str, optional): input delimiter. Defaults to "\t".
        filename_merge (bool, optional): When concatenating,
        create a multi-index with filenames?. Defaults to False.

    Returns:
        pd.Dataframe: concatenated dataframe
    """
    dfs = list()
    input_files = glob.glob(target_glob)
    for fname in sorted(input_files):
        dfs.append(pd.read_csv(fname, sep=input_sep, **kwargs))
    if filename_merge:
        return pd.concat(
            dfs, keys=[Path(x).stem for x in input_files])
    else:
        return pd.concat(dfs, ignore_index=True)


def fig_to_base64(fig):
    """Convert matplot lib fig to code."""
    img = io.BytesIO()
    fig.savefig(
        img, format='png',
        bbox_inches='tight')
    img.seek(0)
    return base64.b64encode(img.getvalue())


def get_circular_stats(input_df: pd.DataFrame, circular_col_name="circ."):
    """Parse Flye output for stats on circularisation.

    Args:
        input_df (pd.DataFrame): Flye stats dataframe
        circular_col_name (str, optional): Col name needed.
        Defaults to "circ.".
    """
    Yes_circ = input_df[input_df[circular_col_name] == "Y"][
        circular_col_name].groupby(level=0).count()
    No_circ = input_df[input_df[circular_col_name] == "N"][
        circular_col_name].groupby(level=0).count()
    merged = pd.concat(
        [Yes_circ, No_circ], axis=1, keys=['Yes_circ', 'No_circ'])
    merged = merged.fillna(0)
    merged = merged.astype(int)
    return merged['Yes_circ']


def run_qc_stats(
        read_stats_glob="stats/*.stats",
        quast_stats_path="assembly_QC.txt",
        flye_stats_glob="flye_stats/*stats.tsv"):
    """Collate data from stats files."""
    # Read stats
    Read_stats = collate_stats(read_stats_glob, input_sep="\t")
    Read_stats_filtered = Read_stats.loc[:, [
        'sample_name', 'read_length', 'mean_quality']]
    # Out table
    Read_median = Read_stats_filtered.groupby(
        'sample_name')['read_length'].median()
    Read_quality = Read_stats_filtered.groupby(
        'sample_name')['mean_quality'].mean()
    Read_MB_sum = Read_stats_filtered.groupby(
        'sample_name')['read_length'].sum()
    Read_count = Read_stats_filtered.groupby(
        'sample_name')['read_length'].count()
    Read_stats_out = pd.concat(
        [Read_count, Read_median, Read_quality, Read_MB_sum], axis=1)
    Read_stats_out.columns = [
        'Read_count', 'Median_read_length', 'Mean_read_quality',
        'Total_read_bases']

    # Quast stats

    Quast_raw_data = pd.read_csv(quast_stats_path, sep='\t', index_col=0)
    Quast_raw_data.index = Quast_raw_data.index.astype(str).str.replace(
        '.medaka', '', regex=True)
    Quast_keep_cols = [12, 13, 14, 15]
    Quast_filtered_data = Quast_raw_data.iloc[:, Quast_keep_cols]
    # Get flye stats

    Flye_stats = collate_stats(flye_stats_glob, filename_merge=True)
    # Remove extension from sample name
    new_index_names = [
        x.replace("_flye_stats", "") for x in
        Flye_stats.index.levels[0].format()]

    Flye_stats.index = Flye_stats.index.set_levels(new_index_names, level=0)

    Flye_cov_mean = Flye_stats[Flye_stats['repeat'] == "N"].groupby(level=0)[
        'cov.'].mean()
    Flye_circ = get_circular_stats(Flye_stats)

    Flye_out = pd.concat([Flye_cov_mean, Flye_circ], axis=1)
    Flye_out.columns = ['Mean_contig_coverage', '#_circular_contigs']

    # merge quast and read stats
    Read_and_assembly = pd.concat(
        [Read_stats_out, Quast_filtered_data, Flye_out], axis=1)
    Read_and_assembly.to_csv('Read_and_assembly_stats.tsv', sep='\t')
    return Read_and_assembly


def gene_plot(gbk_file, **kwargs):
    """Create gene feature plot."""
    color_map = {
        "rep_origin": "yellow",
        "CDS": Colors.cerulean,
        "regulatory": "red",
        "rRNA": Colors.light_cornflower_blue,
        "misc_feature": "lightblue",
    }
    translator = BiopythonTranslator(
        features_filters=(lambda f: f.type not in ["gene", "source"],),
        features_properties=lambda f: {
            "color": color_map.get(f.type, "white")})
    record = translator.translate_record(gbk_file)
    ax, _ = record.plot(
        figure_width=300,
        strand_in_label_threshold=30, **kwargs)
    encoded = fig_to_base64(ax.figure)
    plot = '<pre><img src="data:image/png;base64, {}"></pre>'.format(
        encoded.decode('utf-8'))
    return plot


def gather_sample_files(sample_names):
    """Check files exist for the report."""
    sample_names.sort()
    sample_files = {}
    for sample_name in sample_names:
        variants = os.path.join(
            'variants', sample_name + '.variants.stats')
        depth = os.path.join(
            'total_depth', sample_name + '.total.regions.bed.gz')
        stats = os.path.join(
            'stats', sample_name + '.stats')
        prokka = os.path.join('prokka', sample_name + '.prokka.gbk')
        fwd = os.path.join('fwd', sample_name + '.fwd.regions.bed.gz')
        rev = os.path.join('rev', sample_name + '.rev.regions.bed.gz')
        expected_files = {
            'depth_file': depth,
            'variants_file': variants,
            'prokka': prokka,
            'stats': stats,
            'fwd': fwd,
            'rev': rev}
        final_files = {
            'depth_file': depth,
            'variants_file': variants,
            'prokka': prokka,
            'stats': stats,
            'fwd': fwd,
            'rev': rev}
        for name, file in expected_files.items():
            if os.path.exists(file):
                pass
            else:
                final_files[name] = 'None'
                print('Missing {0} required for report for: {1}'.format(
                    name, sample_name))
        sample_files[sample_name] = final_files
    return sample_files


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser(
        'Visualise mapula output',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--denovo", action="store_true",
        help="Analysis performed de-novo assembly (or variant calling).")
    parser.add_argument(
        "--prokka", action="store_true",
        help="Prokka analysis was performed.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument("--output", help="Report output filename")
    parser.add_argument("--stats", help="directory containing fastcat stats")

    parser.add_argument("--sample_ids", nargs="+")

    args = parser.parse_args()
    report_doc = report.HTMLReport(
        "Bacterial genomes Summary Report",
        ("Results generated through the wf-bacterial-genomes nextflow "
            "workflow provided by Oxford Nanopore Technologies"))
    if not args.denovo:
        report_doc.add_section().markdown(
            "Analysis was completed using an alignment with the provided"
            "reference and medaka was used for variant calling")
    else:
        report_doc.add_section().markdown(
            "As no reference was provided the reads were assembled"
            "and corrected using Flye and Medaka")
        Merged_stats_df = run_qc_stats()
        section = report_doc.add_section()

        section.markdown("## Run summary statistics")
        section.markdown("### Read and assembly statistics")

        section.markdown(
            "This section displays the read and assembly QC"
            " statistics for all the samples in the run.")

        section.table(Merged_stats_df, index=True)

    sample_files = gather_sample_files(args.sample_ids)

    for name, files in sample_files.items():
        section = report_doc.add_section()
        section.markdown("###" + str(name))
        quality_df = pd.read_csv(files['stats'], sep='\t')
        read_qual = fastcat.read_quality_plot(quality_df)
        read_length = fastcat.read_length_plot(quality_df)
        section = report_doc.add_section()
        section.markdown("## Read Quality Control")
        section.markdown(
            "This sections displays basic QC",
            " metrics indicating read data quality.")
        section.plot(
            layout(
                [[read_length, read_qual]],
                sizing_mode="stretch_width")
        )

        section = report_doc.add_section()
        section.markdown("""

        ## Genome coverage
Plots below indicate depth of coverage,
For adequate variant calling depth should be at least 50X in any region.
Forward reads are shown in light-blue, reverse reads are dark grey.
""")
        plots_cover = depthcoverage.depth_coverage(files['depth_file'])
        plots_orient = depthcoverage.depth_coverage_orientation(
            files['fwd'], files['rev'])
        tab1 = Panel(
            child=gridplot(plots_orient, ncols=1),
            title="Coverage traces")
        tab2 = Panel(
            child=gridplot(plots_cover, ncols=1),
            title="Proportions covered")
        cover_panel = Tabs(tabs=[tab1, tab2])
        section.plot(cover_panel)
        if not args.denovo:
            bcfstats.full_report(files['variants_file'], report=section)
        if args.prokka:
            record_dict = SeqIO.to_dict(
                SeqIO.parse(files['prokka'], "genbank"))
            for contig in sorted(record_dict):
                seq_record = record_dict[contig]
                plot = gene_plot(seq_record)
                section = report_doc.add_section()
                section.markdown('##' + str(seq_record.id))
                section.markdown('length:' + str(len(seq_record.seq)))
                section.plot(plot)
    # canned VCF stats report component
    section = report_doc.add_section()
    report_doc.add_section(
        section=scomponents.version_table(args.versions))
    report_doc.add_section(
        section=scomponents.params_table(args.params))
    # Footer section
    section = report_doc.add_section()
    section.markdown('''
### About

**Oxford Nanopore Technologies products are not intended for use for health
assessment or to diagnose, treat, mitigate, cure or prevent any disease or
condition.**

This report was produced using the
[epi2me-labs/wf-bacterial-genomes]
(https://github.com/epi2me-labs/wf-bacterial-genomes). The
workflow can be run using `nextflow epi2me-labs/wf-bacterial-genomes --help`

---
''')

    # write report
    report_doc.write(args.output)


if __name__ == "__main__":
    main()
