# purge_haplotigs

Series of scripts to guide identifying haplotigs in a heterozygous diploid genome assembly (such as using FALCON or FALCON-unzip).

#### The problem

Some parts of a genome may have a very high degree of heterozygosity. This causes contigs for both haplotypes of that part of the genome to be assembled as separate primary contigs, rather than as a contig and an associated haplotig. This can be an issue for downstream analysis.

#### The solution

Identify primary contigs that are haplotigs of other primary contigs, and move them to the haplotig 'pool'. These scripts use mapped read coverage and blast/mummer alignments to guide this process. Dotplots are produced for all flagged contig matches to help the user determine the proper assignment of ambiguous contigs.

## Dependencies

- bash
- bedtools
- Rscript with ggplot2
- perl 
- blast (blastn and makeblastdb)
- make (for mummer)
- gnu parallel (optional but recommended)


## Install

directory structure

```
#!text
-purge_haplotigs
    -bin
    -mummer
    -scripts
```

- pull/clone the git
- compile mummer (this is the standard mummer package but with a modified mummerplot script, the pipeline is setup to run from a local install of mummer)

```
#!bash
cd /path/to/purge_haplotigs/mummer
./configure
make

```

 - add the purge_haplotigs bin to your system PATH (**don't** need to add the mummer directory to the PATH as the scripts will find the relative path to this directory).

```
#!bash
cd /path/to/purge_haplotigs/bin
PATH=$PATH:$PWD

# or add to your .bashrc
# printf "PATH=\$PATH:$PWD\n" >> $HOME/.bashrc
```


## Usage

#### PREPARATION

Map your pacbio subreads, or some decent short reads to your genome assembly (we want a library that produces a nice even coverage; we also want a randombest alignment for multimappers), and then sort the bam. The coverage is used to flag contigs that are likely to be haplotigs (or assembly junk etc).

Index your genome.fasta file with samtools faidx if not already done.

#### STEP 1

Generate a coverage histogram by running the first script. You will also need the bedtools genecov output file for STEP 2.

```
#!text
zz0_coverage_hsitogram.sh aligned.bam
```

#### MANUAL STEP

eyeball the histogram, you should have two peaks, one for haploid level of coverage, the other for diploid level of coverage. Choose cutoffs for low coverage, low point between the peaks, and high coverage.

[Example histogram and choosing cutoffs](https://bitbucket.org/mroachawri/purge_haplotigs/src/cf363f94c00fd865891a0469675d6df4a0813820/examples/example_histogram.png)

#### STEP 2

Run the second script to analyse the coverage on a contig by contig basis--produces a contig `coverage_stats.csv` file with flagged contigs for further analysis.

```
#!text
zz1_analyse_gencov.pl  -i aligned.bam.genecov  -o coverage_stats.csv  -l 30  -m 80  -h 145  [ -j 80  -s 80 ]

REQUIRED:
-i      The output of bedtools genomecov from STEP 1
-o      Output file name (csv format)
-l      The read depth low cutoff (use the histogram to eyeball these cutoffs)
-h      The read depth high cutoff
-m      The low point between the haploid and diploid peaks

OPTIONAL:
-j      Flag contig as "j" (junk) if this percentage or greater of the contig is 
            low/high coverage (default=80)
-s      Flag contig as "s" (suspected haplotig) if this percentage or less of the
            contig is diploid level of coverage (default=80)

```

#### STEP 3

Run the purging pipeline. This script will automatically run a blast search, perform mummer alignments etc. to assess which contigs to reassign and which to keep.

```
#!text
purge_haplotigs.pl  -genome genome.fasta  -coverage coverage_stats.csv

REQUIRED:
-g / -genome        Genome assembly in fasta format. Needs to be indexed with samtools faidx.
-c / -coverage      Contig by contig coverage stats csv file from the previous step.

OPTIONAL:
-t / -threads       Number of worker threads to use. DEFAULT = 4.
-o / -outprefix     Prefix for the curated assembly. DEFAULT = "curated".
-b / -best_match    Percent cutoff for identifying a contig as a haplotig. DEFAULT = 75.
-m / -max_match     Percent cutoff for identifying repetitive contigs. DEFAULT = 250.

```

### ALL DONE! 

You will have five files:

- **<prefix>.fasta**: These are the curated primary contigs
- **<prefix>.haplotigs.fasta**: These are all the haplotigs identified in the initial input assembly. The seq IDs will remain the same but the description field for the contigs will now show the associations, e.g. 
```
#!text
>000000F_003 HAPLOTIG<--000000F_009_HAPLOTIG<--000000F_PRIMARY
```
- **<prefix>.artefacts.fasta**: These are the low/high coverage contigs (identified in STEP 2), and the repetitive contigs. 
- **<prefix>.reassignments.tsv**: These are all the reassignments that were made.
- **<prefix>.contig_associations.log**: This shows the contig association "paths", derived from the contig reassignments. e.g 
```
#!text
000000F,PRIMARY -> 000000F_005,HAPLOTIG
                -> 000000F_009,HAPLOTIG -> 000000F_003,HAPLOTIG
                -> 000000F_010,HAPLOTIG
```

You will also get two directories of dotplots:

- **dotplots_unassigned_contigs:** These are the dotplots that remain unassigned (usually because either their matching contigs was itself reassigned or because it was just under the matching coverage threashold for reassignment).
- **dotplots_reassigned_contigs:** These are the dotplots for the contigs that were reassigned. 

### FURTHER CURATION

You can go through the dotplots and check the assignments. You might wish to move reassigned contigs back into the primary contig pool or purge contigs that were too ambiguous for automatic reassignment. Below are some examples of what might occur. 

[Example: x-axis contig is a haplotig](https://bitbucket.org/mroachawri/purge_haplotigs/src/cf363f94c00fd865891a0469675d6df4a0813820/examples/example_haplotig.png)

[Example: x-axis contig is partially a haplotig](https://bitbucket.org/mroachawri/purge_haplotigs/src/cf363f94c00fd865891a0469675d6df4a0813820/examples/example_partial_haplotig.png)

[Example: x-axis contig is collapsed repeat/assembly junk](https://bitbucket.org/mroachawri/purge_haplotigs/src/cf363f94c00fd865891a0469675d6df4a0813820/examples/example_repetitive_junk_contig.png)