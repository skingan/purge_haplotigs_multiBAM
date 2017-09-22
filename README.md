# Purge Haplotigs

Pipeline to help with curating heterozygous diploid genome assemblies (for instance when assembling using FALCON or FALCON-unzip). 

## Contents

 - [Introduction](#markdown-header-introduction)
 - [Dependencies](#markdown-header-dependencies)
 - [Installation](#markdown-header-installation)
 - [Running Purge Haplotigs](#markdown-header-running-purge-haplotigs)
 - [Further Curation](#markdown-header-further-curation)

## Introduction

#### The problem

Some parts of a genome may have a very high degree of heterozygosity. This causes contigs for both haplotypes of that part of the genome to be assembled as separate primary contigs, rather than as a contig and an associated haplotig. This can be an issue for downstream analysis whether you're working on the haploid or phased-diploid assembly.

#### The solution

Identify pairs of contigs that are syntenic and move one of them to the haplotig 'pool'. The pipeline uses mapped read coverage and blast/lastz alignments to determine which contigs to keep for the haploid assembly. Dotplots are produced for all flagged contig matches, juxtaposed with read-coverage, to help the user determine the proper assignment of any remaining ambiguous contigs. The pipeline will run on either a haploid assembly (i.e. FALCON or FALCON-Unzip primary contigs) or on a phased-diploid assembly (i.e. FALCON-Unzip primary contigs + haplotigs). Here are [two examples](https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Examples) of how Purge Haplotigs can improve a haploid and diploid assembly.

## Dependencies

- Bash
- BEDTools and SAMTools
- blast (blastn and makeblastdb)
- lastz - [website](http://www.bx.psu.edu/~rsharris/lastz/), [github](https://github.com/lastz/lastz)
- Perl (with FindBin, Getopt::Long, Time::Piece, threads, Thread::Semaphore)
- Rscript (with ggplot2 and scales)
- GNU Parallel (optional but highly recommended)


## Installation

Currently only tested on Ubuntu, there is a [Detailed installation](https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Install) example for Ubuntu 16.04 LTS in the wiki.

- Install dependencies, make sure they're in the system PATH
- Pull/clone this git

Either:

- Softlink `purge_haplotigs` to a directory in your system PATH

Or:

- Add the Purge Haplotigs bin directory to your system PATH

```
#!bash
# OPTION 1, SYMLINK 
# to install in $HOME/bin
ln -s /path/to/purge_haplotigs/bin/purge_haplotigs $HOME/bin/purge_haplotigs


# OPTION 2, ADD TO THE SYSTEM PATH
# navigate to the Purge Haploitgs bin directory
cd /path/to/purge_haplotigs/bin

# add it to the system PATH
export PATH=$PATH:$PWD

# optionally add a line to your .bashrc
printf "\nexport PATH=\$PATH:$PWD\n" >> $HOME/.bashrc
```

- That's it! check that it's running ok
```
#!bash
$ purge_haplotigs

USAGE:
purge_haplotigs  readhist,contigcov,purge  [script-specific options]

readhist        First step: generate a read-depth histogram for the genome
contigcov       Second step: get contig-by-contig stats and flag suspect contigs
purge           Third step: run the purge_haplotigs pipeline
```

## Running Purge Haplotigs

There is a [tutorial](https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Tutorial) in the wiki that you can follow which goes into a bit more detail for each step.

#### Preparation

Map your PacBio subreads, or some decent long reads (or even short reads) to your haploid or diploid genome assembly. You'll want to map a library that produces an even coverage and use a 'randombest' alignment for multimappers. Sort and index the bam with `samtools index`. Index your genome.fasta file with `samtools faidx` if you haven't already done so.

#### STEP 1

Generate a coverage histogram by running the first script. This script will produce a histogram png image file for you to look at and a BEDTools 'genomecov' output file that you'll need for STEP 2.

```
#!text
purge_haplotigs  readhist  <bam>

REQUIRED:
<bam>   Sorted bam file of reads aligned to your genome assembly

```

#### MANUAL STEP

You should have a bimodal histogram--one peak for haploid level of coverage, one peak for diploid level of coverage. NOTE: If you're using the phased assembly the diploid peak may be very small. Choose cutoffs for low coverage, low point between the two peaks, and high coverage. Example histograms for choosing cutoffs:

[PacBio subreads on Diploid-phased assembly (Primary + Haplotigs)](examples/phased_coverage_histogram.png)

[Illumina PE reads on Haploid assembly (Primary contigs)](examples/coverage_histogram.png)

#### STEP 2

Run the second script using the cutoffs from the previous step to analyse the coverage on a contig by contig basis. This script produces a contig coverage stats csv file with suspect contigs flagged for further analysis or removal.

```
#!text
purge_haplotigs  contigcov  -i aligned.bam.genecov  -o coverage_stats.csv  -l 30  -m 80  -h 145  [ -j 80  -s 80 ]

REQUIRED:
-i      The bedtools genomecov output that was produced from 'purge_haplotigs readhist'
-o      Choose an output file name (csv format)
-l      The read depth low cutoff (use the histogram to eyeball these cutoffs)
-h      The read depth high cutoff
-m      The low point between the haploid and diploid peaks

OPTIONAL:
-j      Auto-assign contig as "j" (junk) if this percentage or greater of the contig is
        low/high coverage (default = 80, > 100 = don't junk anything)
-s      Auto-assign contig as "s" (suspected haplotig) if this percentage or less of the
        contig is diploid level of coverage (default = 80)

```

#### STEP 3

Run the purging pipeline. This script will automatically run a BEDTools windowed coverage analysis, a blast search, and perform lastz alignments etc. to assess which contigs to reassign and which to keep. The pipeline will make several iterations of purging.

```
#!text
purge_haplotigs  purge  -g genome.fasta  -c coverage_stats.csv -b aligned.sorted.bam

REQUIRED:
-g / -genome        Genome assembly in fasta format. Needs to be indexed with samtools faidx.
-c / -coverage      Contig by contig coverage stats csv file from the previous step.
-b / -bam           Samtools-indexed bam file of aligned reads/subreads to the reference (same
                    one used for generating the read-depth histogram).

OPTIONAL:
-t / -threads       Number of worker threads to use. DEFAULT = 4.
-o / -outprefix     Prefix for the curated assembly. DEFAULT = "curated".
-a / -align_cov     Percent cutoff for identifying a contig as a haplotig. DEFAULT = 70.
-m / -max_match     Percent cutoff for identifying repetitive contigs. DEFAULT = 250.
-wind_len           Length of genome window (in bp) for the coverage track in the dotplots. DEFAULT = 9000.
-wind_step          Step distance for genome windows for coverage track. DEFAULT = 3000.

```

#### All done! 

You will have five files

- **<prefix>.fasta**: These are the curated primary contigs
- **<prefix>.haplotigs.fasta**: These are all the haplotigs identified in the initial input assembly. The seq IDs will remain the same but the description field for the contigs will now show the contig "associations" and purging order, e.g. 
```
#!text
>000000F_003 HAPLOTIG<--000000F_009_HAPLOTIG<--000000F_PRIMARY
```
- **<prefix>.artefacts.fasta**: These are the very low/high coverage contigs (identified in STEP 2). NOTE: you'll probably have mitochondrial/chloroplast/etc. contigs in here with the assembly junk. 
- **<prefix>.reassignments.tsv**: These are all the reassignments that were made, as well as the suspect contigs that weren't reassigned.
- **<prefix>.contig_associations.log**: This shows the contig "associations"/purging order, e.g 
```
#!text
000000F,PRIMARY -> 000000F_005,HAPLOTIG
                -> 000000F_009,HAPLOTIG -> 000000F_003,HAPLOTIG
                -> 000000F_010,HAPLOTIG
```

You will also get two directories of dotplots juxtaposed with read-coverage tracks

- **dotplots_unassigned_contigs:** These are the dotplots for the suspect contigs that remain unassigned.
- **dotplots_reassigned_contigs:** These are the dotplots for the contigs that were reassigned. 

## Further Curation

You can go through the dotplots and check the assignments. You might wish to move reassigned contigs back into the primary contig pool or purge contigs that were too ambiguous for automatic reassignment. Below are some examples of what might occur. 

[A haplotig](examples/haplotig.png)

[Contig is mostly haplotig](examples/haploid_diploid_hemizygous.png) - This example has part of the contig with a diploid level of coverage and part at haploid level with poor alignment to either reference (possibly hemizygous region).

[Haplotig with many tandem repeats](examples/repeat_rich.png)

[Haplotig is palindrome](examples/haplotig_with_palindrome.png)

[Contig circling through string graph 'knots'](examples/repeats_string_graph_short_cut.png) - while this may be a valid string graph path it will likely still confound short read mapping to the haploid assembly.
