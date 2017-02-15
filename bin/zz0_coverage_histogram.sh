#!/bin/bash

# script to drive bedtools genomecov and return a coverage histogram
# modify the bedtools genecov step accordingly if you need to go above 200 depth (i.e. if the diploid peak is > 150)

### SETUP ###

BINDIR=$(dirname `readlink -f $0`)
SCRIPTDIR="$BINDIR/../scripts/"
TMPDIR="tmp_reassign_contigs"

# check dependencies and arguments
EXIT=0

if [[ $# -eq 0 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    printf '\nUsage: zz0_coverage_hsitogram.sh aligned.sorted.bam\n\n'
    exit 1
fi

BEDTOOLS=$(type -p bedtools)
if [[ "$BEDTOOLS" == "" ]]; then
    printf 'ERROR: bedtools not found using "type -p bedtools"\n'
    EXIT=1
fi

RSCRIPT=$(type -p Rscript)
if [[ "$RSCRIPT" == "" ]]; then
    printf 'ERROR: Rscript not found using "type -p Rscript"\n'
    EXIT=1
fi

if ! [[ -f "$SCRIPTDIR/gen_histogram.Rscript" ]]; then
    printf "ERROR: plot script not found at $SCRIPTDIR\n"
    EXIT=1
fi

# die if errors
if [[ $EXIT == 1 ]]; then
    printf "\nOne or more errors encountered, exiting...\n\n"
    exit 1
fi

if ! [[ -d "$TMPDIR" ]]; then
    mkdir $TMPDIR
fi

### RUN ANALYSIS ###

# die if any other failures
set -e

# run genecov
if ! [[ -s "$1.gencov" ]]; then
    printf "bedtools genomecov -ibam $1 -max 200 > $1.gencov\n"
    bedtools genomecov -ibam $1 -max 200 > $1.gencov
else
    printf "SKIP: $1.gencov already exists, skipping bedtools genecov\n"
fi

# prepare the csv
if ! [[ -s "$TMPDIR/$1.histogram.csv" ]]; then
    printf "grep genome $1.gencov | awk '{ print $2 "," $3 }' > $TMPDIR/$1.histogram.csv\n"
    grep genome $1.gencov | awk '{ print $2 "," $3 }' > $TMPDIR/$1.histogram.csv
else
    printf "SKIP: $TMPDIR/$1.histogram.csv already exists, skipping gen histogram .csv step\n"
fi

# make a graph
if ! [[ -s "$1.histogram.png" ]]; then
    printf "$SCRIPTDIR/gen_histogram.Rscript $TMPDIR/$1.histogram.csv $1.histogram.png\n"
    $SCRIPTDIR/gen_histogram.Rscript $TMPDIR/$1.histogram.csv $1.histogram.png
else
    printf "SKIP: $1.histogram.png already exists, skipping gen histogram .png step\n"
fi

# done
printf "\nPipeline finished, the read coverage histogram is saved to $PWD/$1.histogram.png\n"
printf "use this file to determine the depth cutoffs for the next stage\n"


