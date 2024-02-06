#!/usr/bin/env python
"""Accumulate state in checkpoints.json."""
import argparse
import json


def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser(description='Accumulate checkpoint info.')
    parser.add_argument(
        "output_file", help="File to write to")
    parser.add_argument(
        "--checkpoints_file", help="Previous checkpoints file to update", default=None)
    parser.add_argument(
        "--metadata", help="Workflow meta map json")
    parser.add_argument(
        "--checkpoint_data", help="Checkpoint data to process")
    parser.add_argument(
        "--output_definitions", help="Workflow output definitions")

    return parser.parse_args()


def make_checkpoints_file(args):
    """Make the initial checkingtons if this is the 1st time."""
    output_schema = json.load(open(args.output_definitions, 'rb'))

    meta = json.load(open(args.metadata, 'rb'))

    checkpoints = dict(
        files=dict(),
        checkpoints=dict()
    )

    # set up our checkpoints structure to be filled later
    for checkpoint, details in output_schema['checkpoints'].items():
        if details['type'] == 'aggregated':
            checkpoints['checkpoints'][checkpoint] = "incomplete"

        elif details['type'] == 'per-sample':
            samples = dict()

            for sample in meta:
                samples[sample['alias']] = 'incomplete'

            checkpoints['checkpoints'][checkpoint] = samples

    checkpoints['files']['workflow-checkpoints'] = './checkpoints.json'

    return checkpoints


def main():
    """Run the entry point."""
    args = argparser()
    output_schema = json.load(open(args.output_definitions, 'rb'))

    # The 1st checkpoint process doesn't have a file, nextflow passes us null
    # allowing us to create the file
    if args.checkpoints_file:
        checkpoints = json.load(open(args.checkpoints_file, 'rb'))
    else:
        checkpoints = make_checkpoints_file(args)

    # Load the new data for the checkpoint we are dealing with
    new_checkpoint_datas = json.load(open(args.checkpoint_data, 'rb'))

    # Add new data to our existing checkpoints data
    for checkpoint, details in checkpoints['checkpoints'].items():
        for new_checkpoint_data in new_checkpoint_datas:
            if checkpoint == new_checkpoint_data['checkpoint_name']:
                # if it's a sample checkpoint then set the sample in question
                # to complete, else set the overall checkpoint state to complete
                if output_schema['checkpoints'][checkpoint]['type'] == 'per-sample':
                    details[new_checkpoint_data["sample"]] = \
                        new_checkpoint_data['status']
                else:
                    checkpoints['checkpoints'][checkpoint] = \
                        new_checkpoint_data['status']

    # if a sample is not-met for a checkpoint there may be no outputs
    for new_checkpoint_data in new_checkpoint_datas:
        if len(new_checkpoint_data['files']) == 0:
            continue
        for output, path in new_checkpoint_data['files'].items():
            if output_schema['files'][output]['type'] == 'per-sample':

                if output not in checkpoints['files']:
                    checkpoints['files'][output] = dict()

                if new_checkpoint_data['sample'] not in checkpoints['files'][output]:
                    checkpoints['files'][output][new_checkpoint_data['sample']] = path
            else:

                if output not in checkpoints['files']:
                    checkpoints['files'][output] = path

    json.dump(checkpoints, open(args.output_file, 'w'), indent=4)


if __name__ == '__main__':
    main()
