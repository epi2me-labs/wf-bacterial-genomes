"""Create run model."""

import json
import os

from workflow_glue.models.custom import (
    Sample, WorkflowResult)
from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("create_run_model")
    samples = []

    metajson = json.loads(args.metadata)
    for meta in metajson:
        json_file = f"sample_results/{meta['alias']}.json"
        if os.path.exists(json_file):
            with open(json_file, "r") as f:
                samples.append(Sample(**json.loads(f.read())))
        else:
            samples.append(Sample(
                alias=meta["alias"],
                barcode=meta["barcode"],
                sample_type=meta["type"],
                results=dict()
            ))

    workflow = WorkflowResult(
        samples=samples
    )

    workflow.to_json(args.output)

    logger.info(f"Run model written to {args.output}")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("run_model")
    parser.add_argument(
        "--jsons", nargs="+", help="sample results json file(s)."
    )
    parser.add_argument("--metadata")
    parser.add_argument("--output", help="Report output filename")
    return parser
