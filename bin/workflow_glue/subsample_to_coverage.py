#!/usr/bin/env python3

# Script to select the longest reads from a FASTQ file
# to reach a desired genome coverage.

import sys

# Expected genome size in base pairs (e.g. 5 Mb for a bacterial genome)
genome_size = 5_000_000

# Desired coverage depth (e.g. 50×)
target_cov = 50

# Total number of bases to reach the target coverage
bases_needed = genome_size * target_cov

# Container to hold all reads (length, header, sequence, plus line, quality)
reads = []

# Read the FASTQ file in 4-line chunks
with open('filtered.fastq') as f:
    while True:
        hdr = f.readline()
        if not hdr:
            break  # End of file
        seq = f.readline().strip()
        plus = f.readline()
        qual = f.readline()
        reads.append((len(seq), hdr, seq, plus, qual))

# Sort reads by sequence length in descending order (longest reads first)
reads.sort(reverse=True, key=lambda x: x[0])

# Output reads until the accumulated base count meets the coverage target
total = 0
for length, hdr, seq, plus, qual in reads:
    if total >= bases_needed:
        break
    sys.stdout.write(hdr)
    sys.stdout.write(seq + '\n')
    sys.stdout.write(plus)
    sys.stdout.write(qual)
    total += length
