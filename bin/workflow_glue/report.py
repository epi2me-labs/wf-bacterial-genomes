"""Create workflow report."""
import json
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
    parse_mlst
)  # noqa: I100,I202
from .util import get_named_logger, wf_parser  # noqa: ABS101


def collate_flye_stats(dir_path, sample_names, suffix, input_sep="\t", **kwargs):
    """Collate stats files."""
    dfs = []
    valid_samples = []
    for sample_name in sample_names:
        file_path = os.path.join(dir_path, sample_name + suffix)
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep=input_sep, **kwargs)
            dfs.append(df)
            valid_samples.append(sample_name)
    samples_with_missing_files = sorted(set(sample_names) - set(valid_samples))
    if not dfs:
        # not a single sample with a valid assembly
        return None, samples_with_missing_files
    else:
        return (
            pd.concat(
                dfs,
                axis=0,
                keys=valid_samples,
            ),
            samples_with_missing_files,
        )


def get_circular_flye_stats(input_df, circular_col_name="circ."):
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
    flye_stats, samples_with_missing_files = collate_flye_stats(
        flye_dir, sample_names, flye_suffix
    )
    if flye_stats is None:
        return flye_stats, samples_with_missing_files
    flye_stats["circ."] = flye_stats["circ."].str.strip()
    flye_cov_mean = (
        flye_stats[
            flye_stats["repeat"] == "N"].groupby(level=0)["cov."].mean().round(2)
    )
    # flye_cov_mean = flye_cov_mean.reset_index()
    flye_circular = get_circular_flye_stats(flye_stats)
    flye_out = pd.concat([flye_cov_mean, flye_circular["Y"]], axis=1)
    flye_out.columns = ["Mean contig coverage", "# circular contigs"]
    return flye_out, samples_with_missing_files


def gather_sample_files(sample_names, denovo_mode, prokka_mode, isolates_mode, logger):
    """Collect files required for the report per sample and make sure they exist."""
    sample_files = {}
    subdirs_and_suffixes = {
        "total_depth": ["total_depth", "total.regions.bed.gz"],
        "fwd_depth": ["fwd", "fwd.regions.bed.gz"],
        "rev_depth": ["rev", "rev.regions.bed.gz"],
        "variants": ["variants", "variants.stats"],
        "prokka": ["prokka", "prokka.gff"],
        "resfinder": ["resfinder", "resfinder_results.txt"],
        "mlst": ["mlst", "mlst.json"],
        "serotype": ["serotype", "serotype_results.tsv"]
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
            # these may not be produced as errorStrategy is set to "ignore"
            # serotype only produced for salmonella isolates
            elif (
                file_type in ("resfinder", "mlst", "serotype")
                and isolates_mode
            ):
                if not os.path.exists(file):
                    logger.error(
                        f"Isolates file for '{file_type}' missing "
                        f"for sample '{sample_name}' - Check status of analysis."
                        )
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
    "A", "C").

    :param substitution_counts: `pd.DataFrame` with columns `type` and `count`. `type`
        should be of format "X>Y". This is produced as part of the `bcftools stats`
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


def create_report(args, logger):
    """Create and populate Labs report."""
    report = labs.LabsReport(
        "Bacterial Genomes Summary Report",
        "wf-bacterial-genomes",
        args.params,
        args.versions,
    )
    samples = sorted(args.sample_ids)

    if args.stats:
        sample_ids_with_stats = sorted(
            zip(args.sample_ids_with_stats, args.stats), key=lambda x: x[0]
        )
        with report.add_section("Read summary", "Read summary"):
            fastcat.SeqSummary(
                sample_names=tuple([x[0] for x in sample_ids_with_stats]),
                seq_summary=tuple([x[1] for x in sample_ids_with_stats]),
            )

    if not args.denovo:
        html_tags.p(
            "Analysis was completed using an alignment with the provided "
            "reference, and Medaka was used for variant calling."
        )
    else:
        html_tags.p(
            "As no reference was provided the reads were assembled"
            " and corrected using Flye and Medaka."
        )
        stats_table, samples_flye_failed = get_flye_stats(
            sample_names=samples,
            flye_dir="flye_stats",
            flye_suffix="_flye_stats.tsv",
        )
        samples = [
            sample for sample in samples if sample not in samples_flye_failed
        ]
        with report.add_section("Assembly summary statistics", "Assembly"):
            html_tags.p(
                "This section displays the read and assembly QC"
                " statistics for all the samples in the run."
            )
            if stats_table is not None:
                stats_table.index.name = "Sample"
                if samples_flye_failed:
                    html_raw(
                        f"""
                        <b>Info:</b> Flye failed to produce an
                        assembly for the following samples:
                        <b>{", ".join(samples_flye_failed)}</b>.<br>

                        They will be omitted from the rest of the report.
                        """
                    )
                DataTable.from_pandas(stats_table)
            else:
                # not a single sample produced a valid assembly
                html_raw(
                    """
                    <b>Warning:</b> Flye failed to produce an
                    assembly for any of the samples. There are therefore no more results
                    to report.
                    """
                )
                return report

    # Gather stats files for each sample (will be used by the various report sections
    # below)
    sample_files = gather_sample_files(
        samples, args.denovo, args.prokka, args.isolates, logger
        )

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
                for sample in samples:
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
    if args.isolates:
        with report.add_section("Antimicrobial resistance", "AMR"):
            html_raw(
                """
                The contigs were analysed for antimicrobial resistance using <a
                href="https://bitbucket.org/genomicepidemiology/resfinder/src/master/">
                ResFinder</a>. You can use the dropdown menu to select different
                samples.
                """  # noqa
            )
            # add a table with the `resfinder` features for each sample
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for name, files in sample_files.items():
                    with tabs.add_dropdown_tab(name):
                        if files["resfinder"] is None:
                            html_raw(
                                """
                                <b>
                                Resfinder did not work for this sample.
                                Check species and reference parameters are relevant
                                to the sample.</b>
                                """
                            )
                        else:
                            resfinder_df = pd.read_csv(
                                files["resfinder"], sep="\t"
                            ).rename(columns=lambda col: col[0].upper() + col[1:])
                            DataTable.from_pandas(resfinder_df, use_index=False)
        with report.add_section("Multilocus sequence typing", "MLST"):
            html_raw(
                """
                Multilocus sequencing typing was performed on contigs using
                <a href="https://github.com/tseemann/mlst">
                MLST</a>. Typing scheme information is available at
                <a href="https://pubmlst.org/"> PubMLST</a>. You can use the
                dropdown menu to select different samples.
                """
            )
            # Add a table with the `mlst` features for each sample
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for name, files in sample_files.items():
                    with tabs.add_dropdown_tab(name):
                        # use custom func to capitalize the column names because
                        # `str.capitalize()` also changes all-upper-case strings
                        # (e.g. "ID" to "Id")
                        mlst_df = parse_mlst(files["mlst"])
                        if mlst_df is None:
                            html_raw("""
                <b>
                MLST was unable to identify scheme for this sample.
                Please check coverage of sample.</b>
                """)
                        else:
                            mlst_df = mlst_df.rename(
                                columns=lambda col: col[0].upper() + col[1:]
                            )
                            DataTable.from_pandas(mlst_df, use_index=False)
        serotype = False
        # check if we got a serotype for at least one sample
        for files in sample_files.values():
            if files["serotype"]:
                serotype = True
                break
        if serotype:
            with report.add_section("Salmonella serotyping", "Sero."):
                tabs = Tabs()
                with tabs.add_dropdown_menu():
                    for name, files in sample_files.items():
                        with tabs.add_dropdown_tab(name):
                            if files["serotype"] is None:
                                """
                                <b>
                                Serotyping is only available for isolates
                                identified as Salmonella.</b>
                                """
                            else:
                                columns = [
                                    "Predicted serotype",
                                    "Predicted antigenic profile",
                                    "Predicted identification",
                                    "O antigen prediction",
                                    "H1 antigen prediction(fliC)",
                                    "H2 antigen prediction(fljB)",
                                    "Note"
                                ]
                                sero_df = pd.read_csv(
                                    files["serotype"], sep="\t",
                                    usecols=columns
                                )[columns]
                                DataTable.from_pandas(sero_df, use_index=False)

    client_fields = None
    if args.client_fields:
        with open(args.client_fields) as f:
            try:
                client_fields = json.load(f)
            except json.decoder.JSONDecodeError:
                error = "ERROR: Client info is not correctly formatted"

        with report.add_section("Workflow Metadata", "Workflow Metadata"):
            if client_fields:
                df = pd.DataFrame.from_dict(
                    client_fields, orient="index", columns=["Value"])
                df.index.name = "Key"

                # Examples from the client had lists as values so join lists
                # for better display
                df['Value'] = df.Value.apply(
                    lambda x: ', '.join(
                        [str(i) for i in x]) if isinstance(x, list) else x)

                DataTable.from_pandas(df)
            else:
                html_tags.p(error)
    return report


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = create_report(args, logger)
    report.write(args.output)
    logger.info(f"Report written to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("--stats", nargs="*", help="Fastcat per-read stats file(s).")
    parser.add_argument(
        "--sample_ids_with_stats",
        nargs="*",
        help="Sample names in order of per-read stats files.",
    )
    parser.add_argument(
        "--denovo",
        action="store_true",
        help="Analysis performed de-novo assembly (instead of variant calling).",
    )
    parser.add_argument(
        "--prokka", action="store_true", help="Prokka analysis was performed."
    )
    parser.add_argument(
        "--isolates", action="store_true",
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
    parser.add_argument(
        "--client_fields", default=None, required=False,
        help="A JSON file containing useful key/values for display")

    return parser
