#!/usr/bin/env python

import argparse
import glob
import numpy as np
import pandas as pd

from bokeh.layouts import gridplot, layout
from bokeh.models import Panel, Tabs
from bokeh.models.formatters import NumeralTickFormatter
import aplanat
from aplanat import annot, bars, gridplot, hist, lines, points, report, spatial

from parse import parse_bcftools_stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("depth", help="Depth summary file.")
    parser.add_argument("summary", help="Read statistics summary file.")
    parser.add_argument("align_summary", help="Align statistics summary file.")
    parser.add_argument("bcf_stats", help="Output of bcftools stats")
    parser.add_argument("output", help="Report output filename")
    args = parser.parse_args()

    report_doc = report.HTMLReport(
        "Haploid variant calling Summary Report",
        "Results generated through the wf-hap-snp nextflow workflow provided by Oxford Nanopore Technologies")

    report_doc.markdown('''
### Read Quality control
This section displays basic QC metrics indicating read data quality.
''')

    np_blue = '#0084A9'
    np_dark_grey = '#455560'
    np_light_blue = '#90C6E7'

    # read length summary
    seq_summary = pd.read_csv(args.summary, sep='\t')
    total_bases = seq_summary['read_length'].sum()
    mean_length = total_bases / len(seq_summary)
    median_length = np.median(seq_summary['read_length'])
    datas = [seq_summary['read_length']]
    length_hist = hist.histogram(
        datas, colors=[np_blue], bins=100,
        title="Read length distribution.",
        x_axis_label='Read Length / bases',
        y_axis_label='Number of reads',
        xlim=(0, None))
    length_hist = annot.subtitle(
        length_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_length, median_length))

    datas = [seq_summary['acc']]
    mean_q, median_q = np.mean(datas[0]), np.median(datas[0])
    q_hist = hist.histogram(
        datas, colors=[np_blue], bins=100,
        title="Read quality (wrt reference sequence)",
        x_axis_label="Read Quality",
        y_axis_label="Number of reads",
        xlim=(85, 100))
    q_hist = annot.subtitle(
        q_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_q, median_q))

    report_doc.plot(gridplot([[length_hist, q_hist]]))

    report_doc.markdown('''
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
            xs, ys, colors=[np_light_blue, np_dark_grey],
            title="{}: {:.0f}X, {:.1f}% > {}X".format(
                sample, depth.mean(), depth_thresh, depth_lim),
            height=300, width=800,
            x_axis_label='position', y_axis_label='depth',
            ylim=(0,300))
        plots_orient.append(plot)

        # cumulative coverage
        xs = [data['depth'].sort_values(ascending=False)]
        ys = [np.linspace(1, 100, len(data))]
        plot = lines.line(
            xs, ys, colors=[np_light_blue],
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
    report_doc.plot(cover_panel)

    # Footer section
    report_doc.markdown('''
### Variant Calls

This section summarises information from the VCF file found within the workflow output directory.

#### Summary

The following tables and figures are derived from the output of `bcftools stats`.
''')

    vcf_tables = parse_bcftools_stats(args.bcf_stats)
    report_doc.markdown("""
**Variant counts:**
Barcoded samples are represented independently
""")
    df = vcf_tables['SN'].rename(columns={'id':'sample'}) \
        .drop(columns='samples').set_index('sample').transpose()
    report_doc.table(df, index=True)
    report_doc.markdown("**Transitions and tranversions:**")
    df = vcf_tables['TSTV'].rename(columns={'id':'sample'}) \
        .set_index('sample').transpose()
    report_doc.table(df, index=True)
    report_doc.markdown("""
**Substitution types**

Base substitutions aggregated across all samples (symmetrised by pairing)
""")

    sim_sub = {
        'G>A': 'C>T', 'G>C': 'C>G', 'G>T': 'C>A',
        'T>A': 'A>T', 'T>C': 'A>G', 'T>G': 'A>C'}
    def canon_sub(sub):
        b1 = sub[0]
        if b1 not in {'A', 'C'}:
            return canon_sub(sim_sub[sub])
        else:
            return b1, sub[2]

    df = vcf_tables['ST']
    df['canon_sub'] = df['type'].apply(canon_sub)
    df['original'] = df['canon_sub'].apply(lambda x: x[0])
    df['substitution'] = df['canon_sub'].apply(lambda x: x[1])
    df['count'] = df['count'].astype(int)
    df = df[['original', 'substitution', 'count']] \
        .groupby(['original', 'substitution']) \
        .agg(count=pd.NamedAgg(column='count', aggfunc='sum')) \
        .reset_index()

    from bokeh.models import ColorBar, LinearColorMapper
    from bokeh.palettes import Blues9
    from bokeh.plotting import figure
    colors = Blues9[::-1]
    mapper = LinearColorMapper(
        palette=colors, low=min(df['count']), high=max(df['count']))
    p = figure(
        y_range=['C', 'A'], x_range=['A', 'C', 'G', 'T'],
        x_axis_location="above",
        x_axis_label='alternative base',
        y_axis_label='reference base',
        tools="save", toolbar_location='below',
        output_backend="webgl",
        height=225, width=300,
        tooltips=[('sub', '@original>@substitution'), ('count', '@count')])
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.rect(
        source=df, y="original", x="substitution", width=1, height=1,
        fill_color={'field': 'count', 'transform': mapper},
        line_color=None)
    color_bar = ColorBar(
        title='', color_mapper=mapper, label_standoff=10,
        location=(0, 0))
    #p.add_layout(color_bar, 'right')

    report_doc.plot(p)

    report_doc.markdown("""
**Indel lengths**

Insertion and deletion lengths aggregated across all samples.
""")
    df = vcf_tables['IDD']
    df['nlength'] = df['length (deletions negative)'].astype(int)
    df['count'] = df['number of sites'].astype(int)
    # pad just to pull out axes by a minimum
    pad = pd.DataFrame({'nlength':[-10,+10], 'count':[0,0]})
    counts = df.groupby('nlength') \
        .agg(count=pd.NamedAgg(column='count', aggfunc='sum')) \
        .reset_index().append(pad)
    plot = hist.histogram(
        [counts['nlength']], weights=[counts['count']],
        colors = [np_light_blue], binwidth=1,
        title='Insertion and deletion variant lengths',
        x_axis_label='Length / bases (deletions negative)',
        y_axis_label='Count')
    #plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
    report_doc.plot(plot)




    # Footer section
    report_doc.markdown('''
### About

**Oxford Nanopore Technologies products are not intended for use for health assessment
or to diagnose, treat, mitigate, cure or prevent any disease or condition.**

This report was produced using the [epi2me-labs/wf-hap-snp](https://github.com/epi2me-labs/wf-hap-snp).
The workflow can be run using `nextflow epi2me-labs/wf-hap-snp --help`

---
''')

    # write report
    report_doc.write("summary_report.html")

if __name__ == "__main__":
    main()
