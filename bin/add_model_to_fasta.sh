#!/bin/bash
MODEL=$1
FASTA=$2
sed "/>/ s/$/ basecall_model=$MODEL/" $FASTA | bgzip > $FASTA.gz