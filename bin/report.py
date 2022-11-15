#!/usr/bin/env python
"""Create workflow report."""

import argparse
import base64
import io
import os
import re

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


def collate_stats(
        dir_path,
        sample_names,
        suffix,
        input_sep="\t",
        ** kwargs):
    """Collate stats files."""
    dfs = list()
    for sample_name in sample_names:
        file_path = os.path.join(dir_path, sample_name + suffix)
        df = pd.read_csv(file_path, sep=input_sep, **kwargs)
        dfs.append(df)
    return pd.concat(dfs, axis=0, keys=sample_names)


def fig_to_base64(fig):
    """Convert matplot lib fig to code."""
    img = io.BytesIO()
    fig.savefig(
        img, format='png',
        bbox_inches='tight')
    img.seek(0)
    return base64.b64encode(img.getvalue())


def get_circular_stats(input_df, circular_col_name="circ."):
    """Parse flye output for stats on circularisation."""
    circular = input_df.groupby(level=0)[circular_col_name].value_counts()
    circular = circular.unstack()
    circular = circular.fillna(0)

    out_df = pd.DataFrame({
        'N': [0]*len(circular.index),
        'Y': [0]*len(circular.index)},
        index=circular.index)
    for sample in circular.index:
        if 'N' in circular.columns:
            out_df.loc[sample, 'N'] = circular.loc[sample, 'N']
        if 'Y' in circular.columns:
            out_df.loc[sample, 'Y'] = circular.loc[sample, 'Y']
    return out_df


def bp_to_mb(input_num):
    """Convert bases to megabases."""
    return round(input_num / 1000000, 2)


def run_qc_stats(
        sample_names,
        read_stats_dir,
        read_stats_suffix,
        quast_path,
        flye_dir,
        flye_suffix,
):
    """Run QC stats."""
    # read stats
    read_stats = collate_stats(
        read_stats_dir, sample_names, read_stats_suffix)
    read_stats_filtered = read_stats.loc[:, [
        'sample_name', 'read_length', 'mean_quality']]
    # Out table
    read_median = read_stats_filtered.groupby(
        'sample_name')['read_length'].median()
    read_quality = read_stats_filtered.groupby(
        'sample_name')['mean_quality'].mean()
    read_mb_sum = read_stats_filtered.groupby(
        'sample_name')['read_length'].sum()
    read_count = read_stats_filtered.groupby(
        'sample_name')['read_length'].count()
    read_stats_out = pd.concat(
        [read_count,
         read_median,
         read_quality,
         bp_to_mb(read_mb_sum)], axis=1)
    read_stats_out.columns = [
        'Read count', 'Median read length (bp)', 'Mean read quality',
        'Read data (Mb)']

    # quast stats

    quast_raw_data = pd.read_csv(
        os.path.join(quast_path, "transposed_report.tsv"),
        sep='\t',
        index_col=0)
    quast_raw_data.index = quast_raw_data.index.astype(str).str.replace(
        '.medaka', '', regex=True)
    quast_keep_cols = [12, 13, 14, 15]
    quast_filtered_data = quast_raw_data.iloc[:, quast_keep_cols].copy()
    # Get flye stats
    quast_filtered_data.iloc[:, [1, 2, 3]] = bp_to_mb(
        quast_filtered_data.iloc[:, [1, 2, 3]])
    quast_filtered_data.columns = [
        '# contigs',
        'Largest contig (Mb)',
        'Total length (Mb)',
        'Reference length (Mb)']

    # flye stats
    flye_stats = collate_stats(
        flye_dir, sample_names, flye_suffix)
    flye_stats['circ.'] = flye_stats['circ.'].str.strip()
    flye_cov_mean = flye_stats[flye_stats['repeat'] == "N"].groupby(
        level=0)['cov.'].mean()
    # flye_cov_mean = flye_cov_mean.reset_index()
    flye_circular = get_circular_stats(flye_stats)
    flye_out = pd.concat([flye_cov_mean, flye_circular['Y']], axis=1)
    flye_out.columns = ['Mean contig coverage', '# circular contigs']
    # Merge
    merged = pd.merge(
        read_stats_out, quast_filtered_data, left_index=True,
        right_index=True).merge(
            flye_out, left_index=True, right_index=True)
    merged.index.name = None

    return merged


def run_species_stats(species_stats_path, sample_names):
    """Analysis of metaquast data."""
    results = []
    for indexs, sample_name in enumerate(sample_names):
        species_data_path = os.path.join(
            species_stats_path, "blast.res_"+sample_name+"-medaka")
        species_data = pd.read_csv(species_data_path, sep='\t', comment="#")
        top_hit = species_data.iloc[1, :]
        species_split = re.split("\\.", top_hit[1])[2].split(";")
        species = species_split[len(species_split)-1]
        species = re.sub("_", " ", species)
        results.append(
            {'Sample': sample_name,
             'Species': species,
             'Perc_identity': top_hit[2]})
    results_df = pd.DataFrame(results, index=sample_names).iloc[:, 1:3]
    results_df.columns = ['Species ID', 'Identity (%)']
    return results_df


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


def gather_sample_files(sample_names, denovo_mode, prokka_mode):
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
            elif denovo_mode & (name == 'variants_file'):
                pass
            elif prokka_mode & (name == 'prokka'):
                pass
            else:
                final_files[name] = 'None'
                raise FileNotFoundError(
                    'Missing {0} required for report for: {1}'.format(
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
            " and corrected using flye and Medaka")
        merged_stats_df = run_qc_stats(
            sample_names=args.sample_ids,
            read_stats_dir="stats",
            read_stats_suffix=".stats",
            quast_path="quast_stats",
            flye_dir="flye_stats",
            flye_suffix="_flye_stats.tsv")
        section = report_doc.add_section()

        section.markdown("## Run summary statistics")
        section.markdown("### read and assembly statistics")

        section.markdown(
            "This section displays the read and assembly QC"
            " statistics for all the samples in the run.")

        section.table(merged_stats_df, index=True)

        section.markdown("### Species ID")

        section.markdown(
            "This section displays the Species ID as determined by 16S."
            " The table shows the percentage match of the 16S sequence"
            " to the best match of  the SILVA 16S database.")

        species_stats = run_species_stats(
            species_stats_path="quast_stats/quast_downloaded_references",
            sample_names=args.sample_ids)
        section.table(species_stats, index=True)

    sample_files = gather_sample_files(
        args.sample_ids,
        args.denovo,
        args.prokka)

    for name, files in sample_files.items():
        section = report_doc.add_section()
        section.markdown("###" + str(name))
        quality_df = pd.read_csv(files['stats'], sep='\t')
        read_qual = fastcat.read_quality_plot(quality_df)
        read_length = fastcat.read_length_plot(quality_df)
        section = report_doc.add_section()
        section.markdown("## read Quality Control")
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
