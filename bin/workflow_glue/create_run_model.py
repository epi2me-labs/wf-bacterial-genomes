"""Create run model."""

import json
import os

import workflow_glue.results_schema as wf
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
                samples.append(wf.Sample(**json.loads(f.read())))
        else:
            samples.append(wf.Sample(
                alias=meta["alias"],
                barcode=meta["barcode"],
                sample_type=meta["type"],
                results=dict()
            ))

    workflow = wf.WorkflowResult(
        samples=samples
    )

    with open(args.output, 'w') as f:
        f.write(json.dumps(workflow.dict(), indent=4))

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
