#!/usr/bin/env python
"""Create workflow report."""

import argparse

from aplanat import annot, hist, lines, points, report
from aplanat.components import bcfstats
from aplanat.util import Colors
from bokeh.layouts import gridplot
from bokeh.models import Panel, Tabs
import numpy as np
import pandas as pd


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("depth", help="Depth summary file.")
    parser.add_argument("summary", help="Read statistics summary file.")
    parser.add_argument("align_summary", help="Align statistics summary file.")
    parser.add_argument("bcf_stats", help="Output of bcftools stats")
    parser.add_argument("output", help="Report output filename")
    args = parser.parse_args()

    report_doc = report.HTMLReport(
        "Haploid variant calling Summary Report",
        ("Results generated through the wf-hap-snp nextflow "
            "workflow provided by Oxford Nanopore Technologies"))

    section = report_doc.add_section()
    section.markdown("""
### Read Quality control
This section displays basic QC metrics indicating read data quality.
""")

    # read length summary
    seq_summary = pd.read_csv(args.summary, sep='\t')
    total_bases = seq_summary['read_length'].sum()
    mean_length = total_bases / len(seq_summary)
    median_length = np.median(seq_summary['read_length'])
    datas = [seq_summary['read_length']]
    length_hist = hist.histogram(
        datas, colors=[Colors.cerulean], bins=100,
        title="Read length distribution.",
        x_axis_label='Read Length / bases',
        y_axis_label='Number of reads',
        xlim=(0, None))
    length_hist = annot.subtitle(
        length_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_length, median_length))

    # read quality
    datas = [seq_summary['acc']]
    mean_q, median_q = np.mean(datas[0]), np.median(datas[0])
    q_hist = hist.histogram(
        datas, colors=[Colors.cerulean], bins=100,
        title="Read quality (wrt reference sequence)",
        x_axis_label="Read Quality",
        y_axis_label="Number of reads",
        xlim=(85, 100))
    q_hist = annot.subtitle(
        q_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_q, median_q))

    section.plot(gridplot([[length_hist, q_hist]]))

    section = report_doc.add_section()
    section.markdown('''
### Genome coverage
Plots below indicate depth of coverage of the coloured by amplicon pool.
For adequate variant calling depth should be at least 50X in any region.
Forward reads are shown in light-blue, reverse reads are dark grey.
''')
    df = pd.read_csv(args.depth, sep='\t')
    plots_orient = list()
    plots_cover = list()
    depth_lim = 50
    for sample in sorted(df['rname'].unique()):
        bc = df['rname'] == sample
        depth = df[bc].groupby('pos')['depth'].sum()
        depth_thresh = 100*(depth >= depth_lim).sum() / len(depth)

        # fwd/rev
        data = df[bc].groupby('pos').sum().reset_index()  # Is this necessary?
        xs = [data['pos'], data['pos']]
        ys = [data['depth_fwd'], data['depth_rev']]
        plot = points.points(
            xs, ys, colors=[Colors.light_cornflower_blue, Colors.feldgrau],
            title="{}: {:.0f}X, {:.1f}% > {}X".format(
                sample, depth.mean(), depth_thresh, depth_lim),
            height=300, width=800,
            x_axis_label='position', y_axis_label='depth',
            ylim=(0, 300))
        plots_orient.append(plot)

        # cumulative coverage
        xs = [data['depth'].sort_values(ascending=False)]
        ys = [np.linspace(1, 100, len(data))]
        plot = lines.line(
            xs, ys, colors=[Colors.light_cornflower_blue],
            title="{}: {:.0f}X, {:.1f}% > {}X".format(
                sample, depth.mean(), depth_thresh, depth_lim),
            height=300, width=800,
            x_axis_label='Coverage / bases',
            y_axis_label='%age of reference')
        plots_cover.append(plot)

    tab1 = Panel(
        child=gridplot(plots_orient, ncols=1), title="Coverage traces")
    tab2 = Panel(
        child=gridplot(plots_cover, ncols=1), title="Proportions covered")
    cover_panel = Tabs(tabs=[tab1, tab2])
    section.plot(cover_panel)

    # canned VCF stats report component
    section = report_doc.add_section()
    bcfstats.full_report(args.bcf_stats, report=section)

    # Footer section
    section = report_doc.add_section()
    section.markdown('''
### About

**Oxford Nanopore Technologies products are not intended for use for health
assessment or to diagnose, treat, mitigate, cure or prevent any disease or
condition.**

This report was produced using the
[epi2me-labs/wf-hap-snp](https://github.com/epi2me-labs/wf-hap-snp).  The
workflow can be run using `nextflow epi2me-labs/wf-hap-snp --help`

---
''')

    # write report
    report_doc.write(args.output)


if __name__ == "__main__":
    main()
