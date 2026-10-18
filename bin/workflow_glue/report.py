"""Create workflow report."""
from enum import Enum
import json
import os

import dacite
from dominate import tags as html_tags
import ezcharts as ezc  # noqa: I100, I202
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots import util as ezc_util
import pandas as pd
from workflow_glue.models.custom import WorkflowResult
from workflow_glue.report_amr import antibiotic_section
from workflow_glue.report_utils import fill_none

from .parsers import (  # noqa: ABS101
    parse_bcftools_stats_multi
)
from .util import get_named_logger, wf_parser  # noqa: ABS101


def gather_sample_files(sample_names, logger):
    """Collect files required for the report per sample and make sure they exist."""
    sample_files = {}
    subdirs_and_suffixes = {
        "total_depth": ["total_depth", "total.regions.bed.gz"],
        "fwd_depth": ["fwd", "fwd.regions.bed.gz"],
        "rev_depth": ["rev", "rev.regions.bed.gz"],
        "variants": ["variants", "variants.stats"]
    }
    for sample_name in sorted(sample_names):
        files = {}
        for file_type, (subdir, suffix) in subdirs_and_suffixes.items():
            file = os.path.join(subdir, f"{sample_name}.{suffix}")
            if not os.path.exists(file):
                file = None
            files[file_type] = file
        # Handle per_read_stats separately since it's in the sample_data directory
        stats_file = os.path.join(
            "per_read_stats", sample_name, "fastcat_stats", "per-read-stats.tsv.gz")
        if os.path.exists(stats_file):
            files["per_read_stats"] = stats_file
        else:
            files["per_read_stats"] = None
        sample_files[sample_name] = files

    return sample_files


def get_depth_plots(total, fwd, rev, max_contigs):
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

    # plot only N largest contigs
    if max_contigs is not None:
        contig_lengths = depth_df.groupby('ref').agg({
            'end': 'max'
        }).reset_index()
        contig_lengths.columns = ['ref', 'length']
        largest_contigs = contig_lengths.nlargest(max_contigs, 'length')['ref'].tolist()
        depth_df = depth_df[depth_df['ref'].isin(largest_contigs)]

    # make one plot for each ref
    plots = []
    for ref, sub_df in depth_df.groupby("ref"):
        p = ezc.lineplot(
            data=sub_df.eval("mean_pos = (start + end) / 2"),
            x="mean_pos",
            y="depth",
            hue="hue",
            title=ref,
            marker=False
        )
        p._fig.xaxis.axis_label = "Position along reference"
        p._fig.yaxis.axis_label = "Sequencing depth / Bases"
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
        data=plt_data["nlength"], x="nlength", discrete=True, weights=plt_data["count"],
        title="Insertion and deletion lengths"
    )
    # labels, formatting
    p._fig.x_range.start = plt_data["nlength"].min() - 1
    p._fig.x_range.end = plt_data["nlength"].max() + 1
    p._fig.xaxis.axis_label = "Length / bases (deletions negative)"
    p._fig.yaxis.axis_label = "Count"
    return p


def _format_plasmid_label(cluster_id):
    """Format plasmid cluster ID label with line breaks for novel clusters."""
    if "novel_" not in cluster_id:
        return html_tags.b(f"Plasmid_{cluster_id}")
    remaining = cluster_id.replace("novel_", "")
    elements = [html_tags.b("Plasmid_novel_"), html_tags.br()]
    if len(remaining) > 16:
        elements.extend([
            html_tags.b(remaining[:16]),
            html_tags.br(),
            html_tags.b(remaining[16:])
        ])
    else:
        elements.append(html_tags.b(remaining))
    return elements


def plasmid_section_accordion(sample, html_id, fill_na="-"):
    """Create accordion-style plasmid section for a single sample."""
    plasmid_clusters = sample.get_plasmid_clusters(fill_na)
    if not plasmid_clusters:
        return html_tags.b(
            "MOB-suite did not identify plasmid contigs in this assembly.")

    summary = sample.get_plasmid_summary()
    html_tags.p(
        f"MOB-suite identified {summary['plasmid_count']} "
        f"plasmid contig{'s' if summary['plasmid_count'] != 1 else ''} "
        f"out of {summary['total_contigs']} "
        f"total contig{'s' if summary['total_contigs'] != 1 else ''} in "
        f"this assembly, representing "
        f"{summary['distinct_clusters']} distinct plasmid"
        f"{'s' if summary['distinct_clusters'] != 1 else ''}."
    )

    _div = html_tags.div(cls="accordion-item")
    row_num = 0
    for cluster_id in sorted(plasmid_clusters.keys()):
        row_num += 1
        group_data = plasmid_clusters[cluster_id]
        contigs = group_data['rows']
        meta = group_data['metadata']
        # Create accordion header
        _head = html_tags.h2(
            id=f"{html_id}_header_{row_num}",
            style=(
                "border: 1px solid rgba(0,0,0,.125); "
                "border-collapse: collapse; padding:0; margin-bottom:0"
            )
        )
        _button = html_tags.button(
            html_tags.span(
                _format_plasmid_label(cluster_id),
                style="flex: 0 0 auto; max-width: 230px;"
            ),
            html_tags.span(
                html_tags.span(
                    f"{meta['num_contigs']} "
                    f"contig{'s' if meta['num_contigs'] > 1 else ''}",
                    cls="badge bg-primary"),
                html_tags.span(
                    f"{meta['total_length']: ,} bp",  # noqa: E203, E231
                    cls="badge bg-primary") if meta['total_length'] else "",
            ),
            cls="accordion-button collapsed",
            type="button",
            data_bs_toggle="collapse",
            data_bs_target=f"#{html_id}_collapse{row_num}",
            aria_expanded="false",
            aria_controls=f"{html_id}_collapse{row_num}",
            style=(
                "display: grid; align-items: center; "
                "grid-template-columns: 200px 1fr max-content max-content; "
                "grid-gap: 25px"
            )
        )
        _head.add(_button)
        _div.add(_head)

        _div1 = html_tags.div(
            id=f"{html_id}_collapse{row_num}",
            cls="accordion-collapse collapse",
            aria_labelledby=f"{html_id}_header_{row_num}",
            data_bs_parent=f"#{html_id}"
        )
        _div2 = html_tags.div(cls="accordion-body")
        _summary = html_tags.div(style="margin-bottom: 15px;")
        _summary.add(
            html_tags.strong("Replicon type(s): "),
            html_tags.span(
                ", ".join(sorted(meta['rep_types'])) if meta['rep_types']
                else fill_na),
            html_tags.br()
        )
        _summary.add(
            html_tags.strong("Predicted mobility: "),
            html_tags.span(
                ", ".join(sorted(meta['mob_types'])) if meta['mob_types']
                else fill_na),
            html_tags.br()
        )
        _summary.add(
            html_tags.strong("Best database match (Mash distance): "),
            html_tags.span(
                ", ".join(sorted(meta['mash_neighbors'])) if meta['mash_neighbors']
                else fill_na),
            html_tags.br()
        )
        _summary.add(
            html_tags.strong("Predicted host range: "),
            html_tags.span(
                ", ".join(meta['host_ranges']) if meta['host_ranges']
                else fill_na),
            html_tags.br()
        )
        _summary.add(
            html_tags.strong("Observed host range: "),
            html_tags.span(
                ", ".join(meta['observed_host_ranges']) if meta['observed_host_ranges']
                else fill_na)
        )
        _div2.add(_summary)

        # Create table with contig details
        _table = html_tags.table(cls="table table-striped")
        _thead = html_tags.thead()
        _thead.add(
            html_tags.tr(
                html_tags.th("Contig"),
                html_tags.th("Length (bp)"),
                html_tags.th("Replicon type"),
                html_tags.th("Relaxase type"),
                html_tags.th("MPF type"),
                html_tags.th("oriT type")
            )
        )
        _table.add(_thead)
        for contig in contigs:
            _tr = html_tags.tr()
            _tr.add(
                html_tags.td(contig['Contig']),
                html_tags.td(contig['Length (bp)']),
                html_tags.td(contig['Replicon type']),
                html_tags.td(contig['Relaxase type']),
                html_tags.td(contig['MPF type']),
                html_tags.td(contig['oriT type'])
            )
            _table.add(_tr)
        _div2.add(_table)
        _div1.add(_div2)
        _div.add(_div1)
    html_tags.p(
        html_tags.small("Mash distance approaching 0 indicates close identity."),
        style="margin-top: 10px;"
    )
    return _div


def create_report(args, logger):
    """Create and populate Labs report."""
    report = labs.LabsReport(
        "Bacterial Genomes Summary Report",
        "wf-bacterial-genomes",
        args.params,
        args.versions,
        args.wf_version
    )
    samples = sorted(args.sample_ids)

    with open(args.params) as f:
        params_data = json.load(f)

    sample_files = gather_sample_files(
        samples, logger
        )

    dacite_config = dacite.Config(
        cast=[Enum]
    )
    fill_na = "-"

    with report.add_section("Sample summary", "Sample summary"):
        with open(args.results) as f:
            results_dict = json.load(f)
            results = dacite.from_dict(
                data_class=WorkflowResult,
                data=results_dict,
                config=dacite_config
            )
        sample_summaries = []
        for sample in results.samples:
            summary_data = sample.get_summary_data(fill_na)
            sample_summaries.append(summary_data)

        headers = [
            "Barcode",
            "Alias",
            "Read count",
            "Median read length",
            "Contig count",
            "MLST",
            "Species call (ANI | Fraction query)"
        ]

        with html_tags.table(cls="table"):
            with html_tags.thead():
                for h in headers:
                    html_tags.th(h)

            with html_tags.tbody():
                sorted_samples = sorted(
                    sample_summaries,
                    key=lambda x: x['barcode'] or x['alias']
                )

                for sample_summary in sorted_samples:
                    ani_str = (
                        f"{sample_summary['ani']: .3f}"
                        if sample_summary['ani'] is not None else fill_na)
                    f_query_str = (
                        f"{sample_summary['f_query']: .3f}"
                        if sample_summary['f_query'] is not None else fill_na)

                    species_display = (
                        f"{sample_summary['species_call']} ({ani_str} | {f_query_str})"
                        if sample_summary['species_call'] != fill_na
                        else fill_none(sample_summary['species_call'], fill_na))

                    with html_tags.tr():
                        html_tags.td(sample_summary['barcode'])
                        html_tags.td(sample_summary['alias'])
                        html_tags.td(str(sample_summary['n_seqs']))
                        html_tags.td(str(sample_summary['median_len']))
                        html_tags.td(str(sample_summary['contig_count']))
                        html_tags.td(sample_summary['mlst'])
                        html_tags.td(species_display)

        html_tags.small(
            "MLST: Multi-locus sequence typing",
            html_tags.br(),
            "ANI: Average Nucleotide Identity",
            html_tags.br(),
            "Fraction query: Fraction of query genome matched to reference",
        )

    samples_with_stats = []
    for sample_name in samples:
        stats_file = sample_files[sample_name].get("per_read_stats")
        if stats_file:
            samples_with_stats.append((sample_name, stats_file))

    if samples_with_stats:
        with report.add_section("Read summary", "Read summary"):
            html_tags.p(
                "Reads were filtered and summary statistics were generated using ",
                html_tags.a("fastcat", href="https://github.com/epi2me-labs/fastcat"),
                ".",
                html_tags.br(),
                "Use the dropdown menu to view different samples."
            )
            fastcat.SeqSummary(
                sample_names=tuple([x[0] for x in samples_with_stats]),
                seq_summary=tuple([x[1] for x in samples_with_stats]),
            )

    if args.denovo:
        with report.add_section("Assembly summary statistics", "Assembly"):
            html_tags.p(
                "Assembly QC statistics reported by ",
                html_tags.a("Flye", href="https://github.com/mikolmogorov/Flye"),
                " following de novo assembly.",
                html_tags.br(),
                "You can use the dropdown menu to view different samples."
            )
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for sample in results.samples:
                    with tabs.add_dropdown_tab(sample.alias):
                        if sample.has_assembly():
                            flye_display = []
                            for stat in sample.results.assembly.flye_stats:
                                flye_display.append({
                                    'Contig': fill_none(stat.seq_name, fill_na),
                                    'Length (bp)': fill_none(stat.length, fill_na),
                                    'Coverage': fill_none(stat.coverage, fill_na),
                                    'Circular': fill_none(stat.circular, fill_na),
                                    'Repeat': fill_none(stat.repeat, fill_na),
                                    'Multiplicity': fill_none(
                                        stat.multiplicity, fill_na)
                                })
                            df = pd.DataFrame(flye_display)
                            DataTable.from_pandas(df, use_index=False)
                        else:
                            html_tags.b(
                                f"Flye failed to produce an assembly for "
                                f"{sample.alias}."
                            )
    else:  # Reference based
        with report.add_section("Variant calling", "Variants"):
            html_tags.p(
                "The following tables and figures are derived from the output of ",
                html_tags.a(
                    "bcftools stats",
                    href="https://samtools.github.io/bcftools/bcftools.html#stats"
                ),
                "."
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
            html_tags.br()
            html_tags.br()
            html_tags.b("Variant counts:")
            DataTable.from_pandas(
                bcf_stats["SN"].drop(columns="samples").set_index("sample")
            )

            html_tags.br()
            html_tags.b("Transitions and transversions:")
            # transition / transversion table
            DataTable.from_pandas(bcf_stats["TSTV"].set_index("sample"))

            # now the plots (substitution heat map and indel length histogram)
            html_tags.br()
            html_tags.b("Base substitution types and indel lengths:")
            html_tags.br()
            html_tags.p(
                "Base substitutions were symmetrised by pairing.",
                html_tags.br(),
                "You can select different samples from the dropdown menu."
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
                                html_tags.b(
                                    "No substitutions were found for this sample."
                                )
                            if indel_hist is not None:
                                EZChart(indel_hist, "epi2melabs")
                            else:
                                html_tags.b("No indels were found for this sample.")

    if args.plasmid_id:
        with report.add_section("Plasmid identification", "Plasmids"):
            html_tags.p(
                "Plasmid identification was performed using ",
                html_tags.a("MOB-suite", href="https://github.com/phac-nml/mob-suite"),
                html_tags.br(),
                "Use the dropdown menu to view different samples.",
                html_tags.br(),
                "Click on each plasmid cluster to view detailed information "
                "about the contigs."
            )
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for sample in results.samples:
                    with tabs.add_dropdown_tab(sample.alias):
                        accordion_id = (
                            f"accordionPlasmids_{sample.alias.replace(' ', '_')}")
                        with html_tags.div(cls="accordion", id=accordion_id):
                            plasmid_section_accordion(sample, accordion_id, fill_na)

    if args.bakta:
        with report.add_section("Annotations", "Annot."):
            if params_data.get("run_funannotate") and not params_data.get("run_bakta"):
                annot_tool = html_tags.a(
                    "funannotate",
                    href="https://github.com/nextgenusfs/funannotate")
            else:
                annot_tool = html_tags.a(
                    "Bakta", href="https://github.com/oschwengers/bakta")
            html_tags.p(
                "The contigs were annotated with ",
                annot_tool,
                ". ",
                html_tags.br(),
                "You can use the dropdown menu to view different samples."
            )
            # add a table with the `bakta` features for each sample
            tabs = Tabs()
            fill_na = "-"
            with tabs.add_dropdown_menu():
                for sample in results.samples:
                    with tabs.add_dropdown_tab(sample.alias):
                        if sample.has_annotations():
                            # Convert list of Annotation objects to DataFrame
                            annot_data = []
                            for annot in sample.results.assembly.annotations:
                                annot_data.append({
                                    'Contig': fill_none(annot.contig, fill_na),
                                    'ID': fill_none(annot.ID, fill_na),
                                    'Start': fill_none(annot.start, fill_na),
                                    'End': fill_none(annot.end, fill_na),
                                    'Strand': fill_none(annot.strand, fill_na),
                                    'Gene': fill_none(annot.gene, fill_na),
                                    'Product': fill_none(annot.product, fill_na),
                                    'EC_number': fill_none(annot.ec_number, fill_na),
                                })
                            df = pd.DataFrame(annot_data)
                            df = df.rename(columns=lambda col: col[0].upper() + col[1:])
                            DataTable.from_pandas(df, use_index=False)
                        else:
                            html_tags.b(
                                "No annotation data available for this sample."
                            )

    if args.isolates:
        with report.add_section("Antimicrobial resistance prediction", "AMR"):
            html_tags.p(
                "The assembly was analysed for presence of antimicrobial resistance"
                " genes and mutations using ",
                html_tags.a(
                    "ResFinder",
                    href="https://bitbucket.org/genomicepidemiology/resfinder/src/master/"  # noqa: E501
                ),
                ".",
                html_tags.br(),
                "You can use the dropdown menu to select different samples."
            )
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for sample in results.samples:
                    with tabs.add_dropdown_tab(sample.alias):
                        if sample.has_amr():
                            antibiotic_section(
                                sample.results.antimicrobial_resistance
                            )
                        else:
                            html_tags.b(
                                "No AMR genes or mutations detected for this sample.")
            if params_data.get("pointfinder_ignore_indels", False):
                html_tags.small(
                    "PointFinder analysis was run with indels ignored.")
                html_tags.br()
            if params_data.get("pointfinder_ignore_stop_codons", False):
                html_tags.small(
                    "PointFinder analysis was run with stop codons ignored.")

        with report.add_section("Multilocus sequence typing", "MLST"):
            html_tags.p(
                "Multilocus sequencing typing was performed using ",
                html_tags.a("MLST", href="https://github.com/tseemann/mlst"),
                ". Typing scheme information is available at ",
                html_tags.a("PubMLST", href="https://pubmlst.org/"),
                ".",
                html_tags.br(),
                "You can use the dropdown menu to select different samples."
            )
            # Add a table with the `mlst` features for each sample
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for sample in results.samples:
                    with tabs.add_dropdown_tab(sample.alias):
                        if sample.has_mlst():
                            mlst_data = sample.results.sequence_typing
                            # Build the display dictionary
                            mlst_display = {
                                'Scheme': mlst_data.detected_species,
                                'Sequence Type': fill_none(
                                    mlst_data.sequence_type, fill_na),
                            }
                            # Add typing schema alleles as columns
                            if mlst_data.typing_schema:
                                for schema in mlst_data.typing_schema:
                                    mlst_display[schema.schema_identifier] = fill_none(
                                        schema.allele_variant, fill_na
                                    )
                            # Create DataFrame and apply column formatting
                            df = pd.DataFrame([mlst_display])
                            df = df.rename(columns=lambda col: col[0].upper() + col[1:])
                            DataTable.from_pandas(df, use_index=False)
                        else:
                            html_tags.b(
                                "MLST was unable to identify scheme for "
                                "this sample. Please check coverage of sample."
                            )

        # check if we got a serotype for at least one sample
        if any(sample.has_serotype() for sample in results.samples):
            with report.add_section("Salmonella serotyping", "Sero."):
                html_tags.p(
                    "Serotyping was performed using ",
                    html_tags.a("SeqSero2", href="https://github.com/denglab/SeqSero2"),
                    ".",
                    html_tags.br(),
                    "You can use the dropdown menu to select different samples."
                )
                tabs = Tabs()
                with tabs.add_dropdown_menu():
                    for sample in results.samples:
                        with tabs.add_dropdown_tab(sample.alias):
                            if sample.has_serotype():
                                sero_data = sample.results.serotyping
                                sero_display = {
                                    'Predicted serotype': fill_none(
                                        sero_data.predicted_serotype, fill_na
                                    ),
                                    'Predicted antigenic profile': fill_none(
                                        sero_data.predicted_antigenic_profile, fill_na
                                    ),
                                    'O antigen prediction': fill_none(
                                        sero_data.o_antigen_prediction, fill_na
                                    ),
                                    'H1 antigen prediction(fliC)': fill_none(
                                        sero_data.h1_antigen_prediction, fill_na
                                    ),
                                    'H2 antigen prediction(fljB)': fill_none(
                                        sero_data.h2_antigen_prediction, fill_na
                                    ),
                                    'Note': fill_none(sero_data.note, fill_na)
                                }
                                df = pd.DataFrame([sero_display])
                                DataTable.from_pandas(df, use_index=False)
                            else:
                                html_tags.b(
                                    "Serotyping is only available for isolates ",
                                    "identified as Salmonella though MLST."
                                )

    with report.add_section("Genome coverage", "Depth"):
        html_tags.p(
            "The plot below illustrates depth of coverage. For adequate variant "
            "calling, depth should be at least 50X in any region.",
            html_tags.br(),
            "Use the dropdown menu to select different samples."
        )
        max_coverage = params_data.get('max_coverage_plots')
        if max_coverage is not None:
            html_tags.p(
                f"Coverage for up to {max_coverage} longest contigs is displayed "
                f"for each sample."
            )
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for name, files in sample_files.items():
                with tabs.add_dropdown_tab(name):
                    if files["total_depth"] \
                        and files["fwd_depth"] \
                            and files["rev_depth"]:
                        depth_plots = get_depth_plots(
                            files["total_depth"],
                            files["fwd_depth"],
                            files["rev_depth"],
                            params_data.get('max_coverage_plots', None)
                        )
                        for plot in depth_plots:
                            EZChart(plot, "epi2melabs")
                    else:
                        html_tags.b(
                            "No coverage data for sample, "
                            "check input data and assembly quality."
                        )

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
    parser.add_argument(
        "--denovo",
        action="store_true",
        help="Analysis performed de-novo assembly (instead of variant calling).",
    )
    parser.add_argument(
        "--bakta", action="store_true", help="Bakta analysis was performed."
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
        "--plasmid_id", action="store_true",
        help="MOB-suite plasmid identification was performed."
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
    parser.add_argument(
        "--wf_version", default='unknown',
        help="version of the executed workflow")
    parser.add_argument("--results", help="results JSON for reporting")

    return parser
