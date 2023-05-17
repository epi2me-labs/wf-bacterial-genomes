"""Create workflow report."""
import os

from dominate import tags as html_tags
from dominate.util import raw as html_raw
import pandas as pd

import ezcharts as ezc  # noqa: I100, I202
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots import AxisLabelFormatter, util as ezc_util

from .parsers import (  # noqa: ABS101
    parse_bcftools_stats_multi,
    parse_prokka_gff,
)  # noqa: I100,I202
from .util import get_named_logger, wf_parser  # noqa: ABS101


def collate_stats(dir_path, sample_names, suffix, input_sep="\t", **kwargs):
    """Collate stats files."""
    dfs = list()
    for sample_name in sample_names:
        file_path = os.path.join(dir_path, sample_name + suffix)
        df = pd.read_csv(file_path, sep=input_sep, **kwargs)
        dfs.append(df)
    return pd.concat(dfs, axis=0, keys=sample_names)


def get_circular_stats(input_df, circular_col_name="circ."):
    """Parse flye output for stats on circularisation."""
    circular = input_df.groupby(level=0)[circular_col_name].value_counts()
    circular = circular.unstack()
    circular = circular.fillna(0)

    out_df = pd.DataFrame(
        {"N": [0] * len(circular.index), "Y": [0] * len(circular.index)},
        index=circular.index,
    )
    for sample in circular.index:
        if "N" in circular.columns:
            out_df.loc[sample, "N"] = circular.loc[sample, "N"]
        if "Y" in circular.columns:
            out_df.loc[sample, "Y"] = circular.loc[sample, "Y"]
    return out_df


def get_flye_stats(sample_names, flye_dir, flye_suffix):
    """Get Flye stats."""
    # flye stats
    flye_stats = collate_stats(flye_dir, sample_names, flye_suffix)
    flye_stats["circ."] = flye_stats["circ."].str.strip()
    flye_cov_mean = (
        flye_stats[flye_stats["repeat"] == "N"].groupby(level=0)["cov."].mean()
    )
    # flye_cov_mean = flye_cov_mean.reset_index()
    flye_circular = get_circular_stats(flye_stats)
    flye_out = pd.concat([flye_cov_mean, flye_circular["Y"]], axis=1)
    flye_out.columns = ["Mean contig coverage", "# circular contigs"]
    return flye_out


def gather_sample_files(sample_names, denovo_mode, prokka_mode):
    """Collect files required for the report per sample and make sure they exist."""
    sample_files = {}
    subdirs_and_suffixes = {
        "total_depth": ["total_depth", "total.regions.bed.gz"],
        "fwd_depth": ["fwd", "fwd.regions.bed.gz"],
        "rev_depth": ["rev", "rev.regions.bed.gz"],
        "variants": ["variants", "variants.stats"],
        "prokka": ["prokka", "prokka.gff"],
        "resfinder": ["resfinder", "resfinder_results.txt"]
    }
    for sample_name in sorted(sample_names):
        files = {}
        # add the paths to the dict corresponding to the current sample
        for file_type, (subdir, suffix) in subdirs_and_suffixes.items():
            file = os.path.join(subdir, f"{sample_name}.{suffix}")
            # the three depth files should always be present, the variants file only
            # when no de-novo assembly was made and the prokka output only in prokka
            # mode
            if (
                file_type in ("total_depth", "fwd_depth", "rev_depth")
                or (file_type == "variants" and not denovo_mode)
                or (file_type == "prokka" and prokka_mode)
            ):
                if not os.path.exists(file):
                    raise ValueError(
                        f"Required file '{file_type}' missing "
                        f"for sample '{sample_name}'."
                    )
            elif (file_type == "resfinder"):
                if not os.path.exists(file):
                    file = None
            else:
                # this covers the cases when files are not needed (e.g. `variants` when
                # doing a de-novo assembly)
                file = None
            files[file_type] = file
        sample_files[sample_name] = files
    return sample_files


def get_depth_plots(total, fwd, rev):
    """Create depth of coverage line plots from `mosdepth` output `.regions.bed` files.

    There will be one plot for each reference. Each plot has three lines (total,
    forward, reverse reads).
    """
    read_csv_kwargs = dict(
        sep="\t", header=None, names=["ref", "start", "end", "depth"]
    )
    # read the depth files
    depth_df = pd.concat(
        (
            pd.read_csv(total, **read_csv_kwargs).assign(hue="total"),
            pd.read_csv(fwd, **read_csv_kwargs).assign(hue="fwd"),
            pd.read_csv(rev, **read_csv_kwargs).assign(hue="rev"),
        )
    )
    # make one plot for each ref
    plots = []
    for ref, sub_df in depth_df.groupby("ref"):
        p = ezc.lineplot(
            data=sub_df.eval("mean_pos = (start + end) / 2"),
            x="mean_pos",
            y="depth",
            hue="hue",
        )
        p.title = dict(text=ref)
        for series in p.series:
            series.showSymbol = False
        p.legend = dict(orient="vertical", right=10, top="center")
        p.xAxis.name = "Position along reference"
        p.yAxis.name = "Sequencing depth / Bases"
        plots.append(p)
    return plots


def get_substitution_heatmap(substitution_counts):
    """Create heatmap illustrating proportions of substitution types.

    The counts are symmetrised by pairing (i.e. the heatmap will only have two rows:
    'A', 'C').

    :param substitution_counts: `pd.DataFrame` with columns `type` and `count`. `type`
        should be of format 'X>Y'. This is produced as part of the `bcftools stats`
        summary provided by `parser.parse_bcftools_stats_multi()` and can be
        looked up with `"ST"`.
    :returns: `ezcharts.plot.Plot` containing the heatmap.
    """
    # adapted from https://github.com/epi2me-labs/aplanat/blob/56934650dc55748b0b38f7e16cc652d767a3c721/aplanat/components/bcfstats.py#L73  # noqa
    sim_sub = {
        "G>A": "C>T",
        "G>C": "C>G",
        "G>T": "C>A",
        "T>A": "A>T",
        "T>C": "A>G",
        "T>G": "A>C",
    }

    def canon_sub(sub):
        b1 = sub[0]
        if b1 not in {"A", "C"}:
            return canon_sub(sim_sub[sub])
        else:
            return b1, sub[2]

    # wrangle the counts
    df = substitution_counts.copy()
    df["canon_sub"] = df["type"].apply(canon_sub)
    df["original"] = df["canon_sub"].apply(lambda x: x[0])
    df["substitution"] = df["canon_sub"].apply(lambda x: x[1])
    df["count"] = df["count"].astype(int)
    df = (
        df[["original", "substitution", "count"]]
        .groupby(["original", "substitution"])
        .agg(count=pd.NamedAgg(column="count", aggfunc="sum"))
        .reset_index()
    )
    df = df.pivot(index="original", columns="substitution", values="count").rename_axis(
        index="Reference base", columns="Alternative base"
    )
    # normalize to percent
    df = (df * 100 / df.sum().sum()).round(1)

    # draw the heatmap
    p = ezc.heatmap(df.T, vmin=0)
    # format x-axis
    p.xAxis.axisLabel.rotate = 0
    p.xAxis.nameLocation = "center"
    p.xAxis.nameGap = 30
    p.xAxis = [p.xAxis, dict(type="category", position="top")]
    # format y-axis
    p.yAxis.nameLocation = "center"
    p.yAxis.nameGap = 30
    # other formatting
    p.series[0].itemStyle.borderColor = "black"
    p.series[0].itemStyle.borderWidth = 0.5
    p.series[0].label.formatter = "{@[2]} %"
    p.title = dict(text="Substitution types")
    # use our blue for the `visualMap` and hide the slider
    p.visualMap[0].inRange["color"] = ["white", ezc_util.choose_palette()[0]]
    p.visualMap[0].show = False
    return p


def get_indel_length_histogram(indel_lengths):
    """Create a histogram of indel lengths.

    :param indel_lengths: `pd.DataFrame` created as part of the `bcftools stats` summary
        produced by `parser.parse_bcftools_stats_multi()` (can be looked up from the
        resulting `dict` with `"IDD"`).
    :returns: `ezcharts.plot.Plot` containing the histogram.
    """
    # adapted from https://github.com/epi2me-labs/aplanat/blob/56934650dc55748b0b38f7e16cc652d767a3c721/aplanat/components/bcfstats.py#L137  # noqa
    df = indel_lengths[["length (deletions negative)", "number of sites"]].astype(int)
    df.columns = ["nlength", "count"]
    # To draw a histogram with seaborn / ezCharts from pre-binned data we need to create
    # a number range spanning the complete range of x-values. The height of the bars
    # will then be determined by the `weight` parameter.
    plt_data = pd.Series(0, index=range(df["nlength"].min(), df["nlength"].max() + 1))
    plt_data[df["nlength"]] = df["count"]
    plt_data = plt_data.reset_index()
    plt_data.columns = ["nlength", "count"]
    # now plot the histogram with one bar at each position in `nlength` and bar heights
    # corresponding to `count`
    p = ezc.histplot(
        data=plt_data["nlength"], x="nlength", discrete=True, weights=plt_data["count"]
    )
    # title, labels, formatting
    p.title = dict(text="Insertion and deletion lengths")
    p.xAxis.min = plt_data["nlength"].min() - 1
    p.xAxis.max = plt_data["nlength"].max() + 1
    p.xAxis.name = "Length / bases (deletions negative)"
    p.yAxis.name = "Count"
    # make sure we don't use scientific notation for x-axis tick labels
    p.xAxis.axisLabel = dict(formatter=AxisLabelFormatter(sci_notation=False))
    return p


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "Bacterial Genomes Summary Report",
        "wf-bacterial-genomes",
        args.params,
        args.versions,
    )

    if args.stats:
        with report.add_section("Read summary", "Read summary"):
            fastcat.SeqSummary(args.stats)

    if not args.denovo:
        html_tags.p(
            "Analysis was completed using an alignment with the provided "
            "reference, and Medaka was used for variant calling."
        )
        stats_table = None
    else:
        html_tags.p(
            "As no reference was provided the reads were assembled"
            " and corrected using Flye and Medaka."
        )
        stats_table = get_flye_stats(
            sample_names=args.sample_ids,
            flye_dir="flye_stats",
            flye_suffix="_flye_stats.tsv",
        )
        stats_table.index.name = "Sample"
        with report.add_section("Assembly summary statistics", "Assembly"):
            html_tags.p(
                "This section displays the read and assembly QC"
                " statistics for all the samples in the run."
            )
            if stats_table is not None:
                DataTable.from_pandas(stats_table)

    # Gather stats files for each sample (will be used by the various report sections
    # below)
    sample_files = gather_sample_files(
        args.sample_ids, args.denovo, args.prokka)

    with report.add_section("Genome coverage", "Depth"):
        html_tags.p(
            """
            The plot below illustrates depth of coverage. For adequate variant
            calling, depth should be at least 50X in any region. There will be
            several plots if there were multiple references. Use the dropdown menu
            to select different samples.
            """
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for name, files in sample_files.items():
                with tabs.add_dropdown_tab(name):
                    depth_plots = get_depth_plots(
                        files["total_depth"], files["fwd_depth"], files["rev_depth"]
                    )
                    for plot in depth_plots:
                        EZChart(plot, "epi2melabs")

    if not args.denovo:
        with report.add_section("Variant calling", "Variants"):
            html_raw(
                """
            The following tables and figures are derived from the output of <a
            href="https://samtools.github.io/bcftools/bcftools.html#stats">bcftools
            stats</a>.
            """
            )
            # we need a list of variant files and a corresponding list of sample names
            # for `parse_bcftools_stats_multi()` --> extract from the `sample_files`
            # dict
            variant_files = []
            sample_names = []
            for sample_name, files in sample_files.items():
                variant_files.append(files["variants"])
                sample_names.append(sample_name)
            bcf_stats = parse_bcftools_stats_multi(variant_files, sample_names)

            # get the variant counts table
            html_raw(
                """
                <br><br> <b>Variant counts:</b>
                """
            )
            DataTable.from_pandas(
                bcf_stats["SN"].drop(columns="samples").set_index("sample")
            )

            html_raw(
                """
                <br> <b>Transitions and tranversions:</b>
                """
            )
            # transition / transversion table
            DataTable.from_pandas(bcf_stats["TSTV"].set_index("sample"))

            # now the plots (substitution heat map and indel length histogram)
            html_raw(
                """
                <br> <b>Base substitution types and indel lengths:</b><br> Base
                substitutions were symmetrised by pairing. You can select different
                samples from the dropdown menu.
                """
            )
            # get the heatmap of substitution types and histogram of indel lengths for
            # each sample
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for sample in sorted(args.sample_ids):
                    with tabs.add_dropdown_tab(sample):
                        # get substitution heatmap first
                        subst_df = bcf_stats["ST"].query("sample == @sample")
                        subst_heatmap = (
                            get_substitution_heatmap(subst_df)
                            if subst_df["count"].astype(int).sum() > 0
                            else None
                        )
                        # now the indel length histogram
                        indel_hist = None
                        if (
                            "IDD" in bcf_stats
                            and not (
                                indel_lengths_df := bcf_stats["IDD"].query(
                                    "sample == @sample"
                                )
                            ).empty
                        ):
                            indel_hist = get_indel_length_histogram(indel_lengths_df)
                        with Grid():
                            if subst_heatmap is not None:
                                EZChart(subst_heatmap, "epi2melabs")
                            else:
                                html_tags.p(
                                    "No substitutions were found for this sample."
                                )
                            if indel_hist is not None:
                                EZChart(indel_hist, "epi2melabs")
                            else:
                                html_tags.p("No indels were found for this sample.")

    if args.prokka:
        with report.add_section("Annotations", "Annot."):
            html_raw(
                """
                The contigs were annotated with <a
                href="https://github.com/tseemann/prokka">Prokka</a>. You can use the
                dropdown menu to select different samples.
                """
            )
            # add a table with the `prokka` features for each sample
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for name, files in sample_files.items():
                    with tabs.add_dropdown_tab(name):
                        # use custom func to capitalize the column names because
                        # `str.capitalize()` also changes all-upper-case strings
                        # (e.g. "ID" to "Id")
                        prokka_df = parse_prokka_gff(files["prokka"]).rename(
                            columns=lambda col: col[0].upper() + col[1:]
                        )
                        DataTable.from_pandas(prokka_df, use_index=False)
    if args.resfinder:
        with report.add_section("Antimicrobial resistance", "AMR"):
            html_raw(
                """
                The contigs were analysed for antimicrobial resistance using
                <a href="https://bitbucket.org/genomicepidemiology/resfinder/src/master/">
                ResFinder</a>. You can use the
                dropdown menu to select different samples.
                """ # noqa
            )
            # add a table with the `resfinder` features for each sample
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for name, files in sample_files.items():
                    with tabs.add_dropdown_tab(name):
                        if files["resfinder"] is None:
                            html_raw("""
                <b>
                Resfinder did not work for this sample.
                Please check that species and reference parameters are relevant
                to the sample.</b>
                """)
                        else:

                            resfinder_df = pd.read_csv(
                                files["resfinder"], sep='\t').rename(
                                columns=lambda col: col[0].upper() + col[1:]
                            )
                            DataTable.from_pandas(resfinder_df, use_index=False)
    report.write(args.output)
    logger.info(f"Report written to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("--stats", nargs="*", help="Fastcat per-read stats file(s).")
    parser.add_argument(
        "--denovo",
        action="store_true",
        help="Analysis performed de-novo assembly (instead of variant calling).",
    )
    parser.add_argument(
        "--prokka", action="store_true", help="Prokka analysis was performed."
    )
    parser.add_argument(
        "--resfinder", action="store_true",
        help="Resfinder antimicrobial resistance analysis was performed."
    )
    parser.add_argument(
        "--versions",
        required=True,
        help="directory containing CSVs containing name,version.",
    )
    parser.add_argument(
        "--params",
        default=None,
        required=True,
        help="A JSON file containing the workflow parameter key/values",
    )
    parser.add_argument("--output", help="Report output filename")
    parser.add_argument("--sample_ids", nargs="+")
    return parser
