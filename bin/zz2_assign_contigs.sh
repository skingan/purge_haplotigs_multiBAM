#!/bin/bash

# script to drive blastn, and a perl script that will assign contigs and generate dotplots

USAGE="\nUsage: zz2_assign_contigs.sh  stats.csv  genome.fasta\n\n\
NOTE: make sure the genome is indexed with samtools faidx (need genome.fasta.fai)\n\
contig names should only contain alphanumerics, underscores, or hyphens (you'll\n\
need to edit this script otherwise)\n\n"

##### SETUP #####

BINDIR=$(dirname `readlink -f $0`)
SCRIPTDIR="$BINDIR/../scripts/"
TMPDIR="tmp_reassign_contigs"

# check dependencies and args
EXIT=0

if [[ $# -eq 0 ]] || [[ "$1" == "-h" ]] ; then
    printf "$USAGE"
    exit 1
fi

if ! [[ -f $1 ]]; then
    printf "ERROR: $1 is not a file?\n\n$USAGE"
    exit 1
fi

# tempfiles
if ! [[ -d "$TMPDIR" ]]; then
    printf "WARN: temp folder not found, assuming it has been cleaned up--remaking\n"
    mkdir $TMPDIR
fi

BLASTN=$(type -p blastn)
if [[ "$BLASTN" == "" ]]; then
    printf "ERROR blastn not found via type -p blastn\n"
    EXIT=1
fi

# DOTPLOT_MAKER: we assume all of MDP_mummer/nucmer/delta-filter/show-coords are installed if MDP_mummer found
if ! [[ -s "$BINDIR/../dotplot_maker/bin/MDP_mummer" ]]; then
    printf "ERROR dotplotmaker exes not found at $BINDIR/../dotplot_maker/bin/\n"
    EXIT=1
fi

# SCRIPTS: we assume all of the gits scripts are present if analyze_blastn_output.pl is found
if ! [[ -s "$SCRIPTDIR/analyze_blastn_output.pl" ]]; then
    printf "ERROR scripts not found at $SCRIPTDIR\n"
    EXIT=1
fi

### others ???

# die if errors
if [[ $EXIT == 1 ]]; then
    printf '\nOne or more errors encountered, exiting...\n\n'
    exit 1
fi

# die if any fails
set -e

##### START ANALYSIS #####

printf "\n###\n\n"

# make a blast db for later
if ! [[ -d "$TMPDIR/blstdb" ]]; then
    mkdir $TMPDIR/blstdb
fi
if ! [[ -s "$TMPDIR/blstdb/$2.nhr" ]]; then
    printf "\n### Making a blastdb of the reference $2 ###\n"
    printf "makeblastdb -in $2 -dbtype nucl -out $TMPDIR/blstdb/$2\n"
    makeblastdb -in $2 -dbtype nucl -out $TMPDIR/blstdb/$2
else
    printf "SKIP: found blastdb files for $2, skipping generation of blastdb\n"
fi


# get list of suspect contigs
if ! [[ -s "$TMPDIR/suspects.list" ]]; then
    printf "\n### Gathering list of suspected haplotigs from $1 ###\n"
    printf "while IFS='' read -r line || [[ -n \"\$line\" ]]; do\n\
    if [[ \"\$line\" =~ ([a-zA-Z0-9_-]+),[sS], ]]; then\n\
        printf \"\${BASH_REMATCH[1]}\\\\n\" >> $TMPDIR/suspects.list\n\
    elif [[ \"\$line\" =~ ([a-zA-Z0-9_-]+),[jJ], ]]; then\n\
        printf \"\${BASH_REMATCH[1]}\\\\n\" >> $TMPDIR/junk.list\n\
    fi\n\
done < $1\n"
    while IFS='' read -r line || [[ -n "$line" ]]; do
        if [[ "$line" =~ ([a-zA-Z0-9_-]+),[sS], ]]; then
            printf "${BASH_REMATCH[1]}\n" >> $TMPDIR/suspects.list
        elif [[ "$line" =~ ([a-zA-Z0-9_-]+),[jJ], ]]; then
            printf "${BASH_REMATCH[1]}\n" >> $TMPDIR/junk.list
        fi
    done < $1
else
    printf "SKIP: found $TMPDIR/suspects.list, skipping gathering list of suspect contigss from $1\n"
fi

# get the suspect contig seqs
if ! [[ -s "$TMPDIR/suspects.fasta" ]]; then
    printf "\n### getting sequence fastas for suspected haplotigs ready for analysis ###\n"
    printf "$SCRIPTDIR/returnseq.pl -f $2 -l $TMPDIR/suspects.list > $TMPDIR/suspects.fasta\n"
    $SCRIPTDIR/returnseq.pl -f $2 -l $TMPDIR/suspects.list > $TMPDIR/suspects.fasta
else
    printf "SKIP: found suspects.fasta, skipping gathering suspect contig seqs\n"
fi

# do the all-v-all blastn
if ! [[ -s "$TMPDIR/suspects.blastn.gz" ]]; then
    printf "\n### running blastn on suspects.fasta vs genome.fasta ###\n"
    printf "blastn -query $TMPDIR/suspects.fasta -db $TMPDIR/blstdb/$2 -outfmt 6 -num_alignments 3 -evalue 0.0000001 -num_threads 4 |  awk ' $1 != $2 && $4 > 500 { print } ' | gzip - > $TMPDIR/suspects.blastn.gz\n"
    blastn -query $TMPDIR/suspects.fasta -db $TMPDIR/blstdb/$2 -outfmt 6 -num_alignments 3 -evalue 0.0000001 -num_threads 4 |  awk ' $1 != $2 && $4 > 500 { print } ' | gzip - > $TMPDIR/suspects.blastn.gz
else
    printf "SKIP: found $TMPDIR/suspects.blastn.gz, skipping blast search of suspect contigs\n"
fi

# analyse the blastn output, gen dotplots for reciprocal best hits
if ! [[ -s suspect_contig_reassign.tsv ]]; then
    printf "\n### running blastn output analysis, generating dotplots and guessing assignment for suspect contigs\n"
    printf "### this step will take a while...\n"
    printf "$SCRIPTDIR/analyze_blastn_output.pl -f $2 -b $TMPDIR/suspects.blastn.gz -t suspect_contig_reassign.tsv 2> $TMPDIR/analyse_blastn.log\n"
    $SCRIPTDIR/analyze_blastn_output.pl -f $2 -b $TMPDIR/suspects.blastn.gz -t suspect_contig_reassign.tsv 2> $TMPDIR/analyse_blastn.log
else
    printf "SKIP: found suspect_contig_reassign.tsv, skipping blastn analysis and dotplot generation\n"
fi

printf "\n###\n\n    Table is saved in the following format (check and update the table, feed into next step of pipeline).\n";
printf "    Reassign_keys: h = haplotig, r = repeat/assembly junk, c = crop, ? = ¯\\_(^_^)_/¯\n";
printf "    Only 'h', 'r', and 'c' flagged contigs will be reassigned.\n\n"
printf "    For cropping: positions are 1-indexed, can use \"start\" and \"end\" for start/end of contig\n";
printf "    NOTE: The crop region you specify here is the region you wish to KEEP as a primary contig\n"
printf "    (the rest of the contig will be output with the haplotigs).\n\n";
printf "#SUSPECT_CONTIG    TOP_MATCH  SECOND_MATCH  MAXMATCHCOV      BESTMATCHCOV      REASSGIN_KEY  CROP_START  CROP_END\n";
printf "contig1            hit1       hit2          10.00            5.00              \n";
printf "contig2            hit1       hit2          95.00            93.00             h\n";
printf "contig3            hit1       hit2          400.00           98.00             r\n";
printf "contig4            hit1       hit2          50.00            50.00             c             50000       end\n\n";







