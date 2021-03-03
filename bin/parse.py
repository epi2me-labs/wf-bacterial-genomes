import re
import sys

import pandas as pd

def split_blocks(fname):
    """Split lines of a file into blocks based on comment lines."""
    comment = list()
    data = list()
    state = None
    with open(fname, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('#'):
                if state != 'comment' and state is not None:
                    yield comment, data
                    comment, data = list(), list()
                state = 'comment'
                comment.append(line.strip('# ').rstrip())
            else:
                state = 'data'
                data.append(line.rstrip())
        yield comment, data


def parse_bcftools_stats(fname):
    tables = dict()
    with open(fname, 'r') as fh:
        for comment, data in split_blocks(fname):
            fields = [x.rstrip() for x in re.split('\[\d+\]', comment[-1])]
            section, fields = fields[0], fields[1:]
            rows = list()
            for d in data:
                items = d.split('\t')
                if items[0] != section:
                    raise ValueError("first data field not equal to section key")
                rows.append(items[1:])
            tables[section] = pd.DataFrame(rows, columns=fields)
    # now some special handling
    SN = tables['SN']
    SN['key'] = SN['key'].apply(lambda x: x.replace('number of ', '').rstrip(':'))
    tables['SN'] = pd.DataFrame(SN.pivot(index='id', columns='key', values='value').to_records())
    return tables


if __name__ == "__main__":
    tables = parse_bcftools_stats(sys.argv[1])
    print(tables['SN'])
