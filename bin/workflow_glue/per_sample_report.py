"""Create a report for each sample."""
import json
import os

from dominate import tags as html_tags
from ezcharts.components.reports import labs
import pandas as pd

from .collect_results import gather_sample_files  # noqa: ABS101
from .parsers import parse_mlst  # noqa: ABS101
from .process_resfinder_iso import (  # noqa: ABS101
    get_acquired_data,
    get_point_data
)
from .util import get_named_logger, wf_parser  # noqa: ABS101


def convert_bp(size):
    """Convert bp values to appropriate metric prefix."""
    size = float(str(size).replace(",", ""))
    for x in ["bp", "Kb", "Mb", "Gb", "Tb"]:
        if size < 1000.0:
            if x == "bp":
                return "{:.0f} {}".format(size, x)
            else:
                return "{:.2f} {}".format(size, x)
        size /= 1000.0
    return size


def lead_section(text):
    """From a dict of key value pairs make a text section."""
    _div = html_tags.div()
    for key, value in text.items():
        if key == 'header':
            _div.add(html_tags.h3(value))
        else:
            _div.add(
                html_tags.small(key.replace("_", " ").upper(), cls="text-muted"),
                html_tags.h5(value if value is not None else "None")
            )
    return _div


def get_run_summary(files, reference=None):
    """Get run summary statistics."""
    run_dict = {}
    total_yield = 0
    median_read_length = 0
    median_read_length = 0
    if files["fastcat"]:
        read_df = pd.read_csv(files["fastcat"], sep="\t")
        total_yield = read_df["read_length"].sum()
        median_read_length = read_df["read_length"].median()
        median_read_quality = read_df["mean_quality"].median()

    taxon = "Unknown"
    if files["mlst"]:
        mlst_df = parse_mlst(files["mlst"])
        if mlst_df is not None:
            taxon = mlst_df.loc[0, "scheme"]

    if reference is None:
        run_dict["taxon"] = taxon.capitalize()
    else:
        run_dict["reference"] = os.path.basename(reference)
    run_dict["total_yield"] = convert_bp(total_yield),
    run_dict["median_read_length"] = convert_bp(median_read_length),
    run_dict["median_read_quality"] = median_read_quality
    return run_dict


def flye_section(flye_file):
    """Extract information on denovo assembly."""
    df = pd.read_csv(flye_file, sep="\t")
    contig_num = df.shape[0]
    total_yield = convert_bp(df["length"].sum())
    circular_num = len(df[df["circ."] == "Y"])

    _div = html_tags.div()
    _table = html_tags.table(cls="table table-striped")
    _thead = html_tags.thead()
    _thead.add(
        html_tags.tr(
            html_tags.th("# Contigs"),
            html_tags.th("Total length"),
            html_tags.th("# Circular contigs"),
        )
    )
    _table.add(_thead)
    _tr = html_tags.tr()
    _tr.add(html_tags.td(contig_num))
    _tr.add(html_tags.td(total_yield))
    _tr.add(html_tags.td(circular_num))
    _table.add(_tr)
    _div.add(_table)
    return _div


def ref_section(total):
    """Return coverage information on reference based assembly."""
    depth_df = pd.read_csv(total, sep="\t", names=["ref", "start", "end", "depth"])
    coverage_dict = dict()
    thresholds = [30, 50]
    for ref, depths in depth_df.groupby("ref"):
        ref_dict = dict()
        for t in thresholds:
            total_length = 0
            above_threshold = 0
            for index, row in depths.iterrows():
                total_length += row["end"] - row["start"]
                if row["depth"] >= t:
                    above_threshold += row["end"] - row["start"]
            above_threshold_pct = above_threshold / total_length * 100
            if "total_len" not in ref_dict:
                ref_dict["total_len"] = convert_bp(total_length)
            ref_dict[t] = round(above_threshold_pct, 3)
        coverage_dict[ref] = ref_dict

    # HTML table
    _div = html_tags.div()
    _table = html_tags.table(cls="table table-striped")
    _thead = html_tags.thead()
    _thead.add(
        html_tags.tr(
            html_tags.th("Contig"),
            html_tags.th("Total length"),
            [html_tags.th(f"% Coverage at {t}x") for t in thresholds]
        )
    )
    _table.add(_thead)
    for k, v in coverage_dict.items():
        _tr = html_tags.tr()
        _tr.add(html_tags.td(k))
        _tr.add(html_tags.td(v["total_len"]))
        _tr.add([html_tags.td(v[t]) for t in thresholds])
        _table.add(_tr)
    _div.add(_table)
    return _div


def amr_section(resfinder_json, html_id):
    """Parse resfinder JSON for accordion style table."""
    with open(resfinder_json) as fh:
        resfinder_data_raw = json.load(fh)
    acquired_data = get_acquired_data(resfinder_data_raw)
    point_data = get_point_data(resfinder_data_raw)
    if (not acquired_data and not point_data):
        return html_tags.p("No AMR genes detected in sample.")
    _div = html_tags.div(cls="accordion-item")
    # Point mutations
    row = 0   # Fault if gene name contains "'" (e.g aac(2')-Ic
    for gene, evidence in point_data.items():
        if not evidence:
            continue
        row += 1
        drugs = {drug.capitalize() for mut in evidence for drug in mut["drugs"]}
        _head = html_tags.h2(id=f"{row}", style="border: 1px solid rgba(0,0,0,.125);\
                            border-collapse: collapse;\
                             padding:0;\
                             margin-bottom:0")
        _button = html_tags.button(
            html_tags.span(html_tags.b(gene)),
            html_tags.span(
                [html_tags.span(d, cls="badge bg-dark me-1") for d in drugs]
            ),
            html_tags.span(
                "PointFinder",
                html_tags.span(
                    len(evidence),
                    cls="position-absolute badge rounded-pill bg-danger",
                    style="top:8px"
                    ),
                cls="badge bg-primary"),
            cls="accordion-button collapsed",
            type="button",
            data_bs_toggle="collapse",
            data_bs_target=f"#collapse{row}",
            aria_expanded="false",
            aria_controls=f"collapse{row}",
            style="display: grid; \
                    align-items: center;\
                    grid-template-columns: 150px 1fr max-content max-content;\
                    grid-gap: 25px"
            )
        _head.add(_button)
        _div.add(_head)
        _div1 = html_tags.div(
            id=f"collapse{row}",
            cls="accordion-collapse collapse",
            fr=f"{row}",
            aria_labelledby=f"{row}",
            data_bs_parent=f"#{html_id}")
        _div2 = html_tags.div(cls="accordion body")
        _table = html_tags.table(cls="table table-striped")
        _thead = html_tags.thead()
        _thead.add(
            html_tags.tr(
                html_tags.th("Antimicrobial resistance"),
                html_tags.th("Amino Acid"),
                html_tags.th("Nucleotide"),
            )
        )
        _table.add(_thead)
        for mutation in evidence:
            _tr = html_tags.tr()
            _tr.add(
                html_tags.td("; ".join(d.capitalize() for d in mutation["drugs"])),
                html_tags.td(mutation["aa"]),
                html_tags.td(mutation["nuc"].upper())
            )
            _table.add(_tr)
        _div2.add(_table)
        _div1.add(_div2)
        _div.add(_div1)

    # Acquired resistance
    for gene, evidence in acquired_data.items():
        row += 1
        _head = html_tags.h2(id="f{row}", style="border: 1px solid rgba(0,0,0,.125);\
                            border-collapse: collapse;\
                            padding:0;\
                            margin-bottom:0")
        _button = html_tags.button(
            html_tags.span(html_tags.b(gene)),
            html_tags.span(
                {html_tags.span(
                    d.capitalize(),
                    cls="badge bg-dark me-1"
                    ) for d in evidence["drugs"]}
                ),
            html_tags.span(
                "ResFinder",
                html_tags.span(
                    1,  # TODO add evidence as list for multiple shits of gene
                    cls="position-absolute badge rounded-pill bg-danger",
                    style="top:8px"
                    ),
                cls="badge bg-secondary"),
            cls="accordion-button collapsed",
            type="button",
            data_bs_toggle="collapse",
            data_bs_target=f"#collapse{row}",
            aria_expanded="false",
            aria_controls=f"collapse{row}",
            style="display: grid;\
                grid-template-columns: 150px 1fr max-content max-content;\
                align-items: center;\
                grid-gap: 25px")
        _head.add(_button)
        _div.add(_head)
        _div1 = html_tags.div(
            id=f"collapse{row}",
            cls="accordion-collapse collapse",
            fr=f"{row}",
            aria_labelledby=f"{row}",
            data_bs_parent=f"#{html_id}"
            )
        _div2 = html_tags.div(cls="accordion body")
        _table = html_tags.table(cls="table table-striped")
        _thead = html_tags.thead()
        _thead.add(
            html_tags.tr(
                html_tags.th("Antimicrobial resistance"),
                html_tags.th("Identity"),
                html_tags.th("Coverage")
            )
        )
        _table.add(_thead)
        _tr = html_tags.tr()
        _tr.add(
            html_tags.td("; ".join(d.capitalize() for d in evidence["drugs"])),
            html_tags.td(evidence["identity"]),
            html_tags.td(evidence["coverage"])
        )
        _table.add(_tr)
        _div2.add(_table)
        _div1.add(_div2)
        _div.add(_div1)
    return _div


def mlst_section(mlst_file):
    """Extract mlst results."""
    mlst_df = parse_mlst(mlst_file)
    if mlst_df is None:
        return html_tags.div(html_tags.p(
            "MLST was unable to identify a typing scheme for this sample. "
            "Please check coverage of sample."
            ))
    col_names = mlst_df.columns
    row_data = mlst_df.iloc[0].values
    # C apitalise all non-gene columns
    col_names = [x.capitalize() if i < 3 else x for i, x in enumerate(col_names)]
    # HTML section
    _div = html_tags.div()
    _table = html_tags.table(cls="table table-striped")
    _thead = html_tags.thead()
    _thead.add(html_tags.tr([html_tags.th(c) for c in col_names]))
    _table.add(_thead)
    _tr = html_tags.tr()
    for cell in row_data:
        _tr.add(html_tags.td(cell))
    _table.add(_tr)
    _div.add(_table)
    return _div


def serotype_section(serotype_file):
    """Extract serotyping results."""
    columns = [
        "Predicted serotype",
        "Predicted antigenic profile",
        "Predicted identification",
        "O antigen prediction",
        "H1 antigen prediction(fliC)",
        "H2 antigen prediction(fljB)",
        "Note"
    ]
    # Check for file is made before section created in main()
    sero_df = pd.read_csv(
        serotype_file, sep="\t",
        usecols=columns
    )[columns]
    # Done manually to allow non-scrollable table for PDF output
    row_data = sero_df.iloc[0].values
    _div = html_tags.div()
    _table = html_tags.table(cls="table table-striped")
    _thead = html_tags.thead()
    _thead.add(html_tags.tr([html_tags.th(c) for c in columns]))
    _table.add(_thead)
    _tr = html_tags.tr()
    for cell in row_data:
        _tr.add(html_tags.td(cell))
    _table.add(_tr)
    _div.add(_table)
    return _div


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        f"{args.sample_alias} | Isolate Sequencing Report",
        "wf-bacterial-genomes",
        args.params,
        args.versions,
    )

    files = gather_sample_files(args.sample_alias, args.data_dir)

    lead_summary = lead_section(dict(
        sample=args.sample_alias,
        barcode=args.sample_barcode,
        run=args.wf_session,
        version=args.wf_version,
    ))
    with open(args.params) as fh:
        params_data = json.load(fh)

    run_summary_dict = get_run_summary(
        files,
        params_data["reference"]
        )

    run_summary = lead_section(run_summary_dict)

    with report.add_section('Workflow details', 'Workflow', True):
        with html_tags.div(cls="row"):
            html_tags.div(lead_summary, cls="col-sm-12", style="float: left; width:50%")
            html_tags.div(run_summary, cls="col-sm-12", style="float: right; width:50%")

    with report.add_section("Multilocus sequencing typing (MLST)", "MLST"):
        with html_tags.div(cls="row"):
            if files["mlst"]:
                mlst_section(files["mlst"])
            else:
                html_tags.p("MLST was unable to be perfomed.")

    if files["serotype"]:
        with report.add_section("Salmonella serotyping", "Serotype"):
            with html_tags.div(cls="row"):
                serotype_section(files["serotype"])

    with report.add_section('Antimicrobial resistance prediction', 'AMR', True):
        with html_tags.div(cls="accordion", id="accordionTable"):
            if files["amr"]:
                html_tags.p(
                    "Analysis was performed using ResFinder. "
                    "Click on each gene for more information."
                )
                amr_section(files["amr"], "accordionTable")
            else:
                html_tags.p("No AMR genes detected for specific species.")

    with report.add_section("Assembly QC", "Assembly", True):
        with html_tags.div(cls="row"):
            if not args.denovo:
                html_tags.p(
                    "Analysis was completed using an alignment with the provided "
                    "reference, and Medaka was used for variant calling."
                )
                if files["depth"]:
                    ref_section(files["depth"])
                else:
                    html_tags.p(
                        "No coverage information, please check input data quality."
                        )
            else:
                html_tags.p(
                    "As no reference was provided the reads were assembled"
                    " and corrected using Flye and Medaka."
                )
                if files["flye"]:
                    flye_section(files["flye"])
                else:
                    html_tags.p(
                        """No denovo assembly for sample.
                         Please check input data quality."""
                        )

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
    parser.add_argument("--sample-alias", required=True)
    parser.add_argument("--sample-barcode", required=True)
    parser.add_argument("--data_dir", required=True, help="Analysis results directory")
    parser.add_argument("--wf-session", required=True)
    parser.add_argument("--wf-version", required=True)
    return parser
