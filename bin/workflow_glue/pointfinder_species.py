#!/usr/bin/env python
"""Obtain pointfinder species option from mlst scheme."""

import json
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run entry point."""
    logger = get_named_logger("pointfinder_species")
    """Extract mlst scheme and assign pointfinder species."""
    pointfinder_dict = {
        "campylobacter_nonjejuni_7": "campylobacter",
        "campylobacter_nonjejuni_8": "campylobacter",
        "campylobacter_nonjejuni": "campylobacter",
        "campylobacter_nonjejuni_6": "campylobacter",
        "campylobacter_nonjejuni_3": "campylobacter",
        "campylobacter_nonjejuni_5": "campylobacter",
        "campylobacter": "campylobacter",
        "campylobacter_nonjejuni_4": "campylobacter",
        "campylobacter_nonjejuni_2": "campylobacter",
        "efaecium": "enterococcus faecium",
        "efaecalis": "enterococcus faecalis",
        "neisseria": "neisseria gonorrhoeae",
        "senterica_achtman_2": "salmonella",
        "ecoli": 'escherichia_coli',
        "klebsiella": "klebsiella",
        "koxytoca": "klebsiella",
        "kaerogenes": "klebsiella",
        "saureus": "staphylococcus aureus",
        "helicobacter": "helicobacter pylori",
        "mycobacteria_2": "mycobacterium tuberculosis",
    }
    with open(args.mlst_json) as f:
        data = json.load(f)
    pointfinder_species = pointfinder_dict.get(data[0]["scheme"], "other")
    logger.info("Pointfinder species identified.")
    sys.stdout.write(pointfinder_species)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("pointfinder_species")
    parser.add_argument(
        "--mlst_json",
        help="MLST json results file")
    return parser
