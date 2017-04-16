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


## Install

directory structure

```
#!text
-purge_haplotigs
    -bin
    -dotplot_maker
        -bin
        -MDP_MUMer
    -scripts
```

- pull/clone the git
- install MDP_mummer (modified mummer package: I made some fixes and changed some defaults to give nicer looking dotplots, otherwise it's essentially just the mummer package)

```
#!bash
cd /path/to/purge_haplotigs/dotplot_maker/MDP_MUMer
./install.sh        # tested on Ubuntu 16-LTS
```

 - add the purge_haplotigs bin to your system PATH (**don't** need to add `dotplot_maker/bin` to the PATH as the scripts will find the relative path to this directory).

```
#!bash
cd /path/to/purge_haplotigs/bin
PATH=$PATH:$PWD
```


## Usage

#### PREP

Map your pacbio subreads, or some decent short reads to your genome assembly (we want a library that produces a nice even coverage; we also want a randombest alignment for multimappers), and then sort the bam. The coverage is used to flag contigs that are likely to be haplotigs (or assembly junk etc).

Index your genome.fasta file with samtools faidx if not already done.

#### STEP 1

Generate a coverage histogram by running the first script. You will also need the bedtools genecov output file for STEP 2.

```
#!text
zz0_coverage_hsitogram.sh aligned.sorted.bam
```

#### MANUAL STEP

eyeball the histogram, you should have two peaks, one for haploid level of coverage, the other for diploid level of coverage. Choose cutoffs for low coverage, low point between the peaks, and high coverage.

[Example histogram and choosing cutoffs](https://bitbucket.org/mroachawri/purge_haplotigs/src/cf363f94c00fd865891a0469675d6df4a0813820/examples/example_histogram.png)

#### STEP 2

Run the second script to analyse the coverage on a contig by contig basis--produces a contig `coverage_stats.csv` file with flagged contigs for further analysis.

```
#!text
zz1_analyse_gencov.pl  -i genecov.out  -o coverage_stats.csv  -l 30  -m 80  -h 145  [ -j 80  -s 80 ]

REQUIRED:
-i      The output of bedtools genomecov from STEP 1
-o      Output file name (csv format)
-l      The read depth low cutoff (use the histogram to eyeball these cutoffs)
-h      The read depth high cutoff
-m      The low point between the haploid and diploid peaks

OPTIONAL:
-j      Flag contig as "j" (junk) if this percentage or greater of the contig is 
            low/high coverage (default=80)
-s      Flag contig as "s" (suspected haplotig) if this percentage or more of the
            contig is haploid level of coverage (default=80)

```

#### STEP 3

Run the purging pipeline. This script will automatically run a number of steps to identify contig matches for 'suspect' contigs, perform mummer alignments, and will attempt to guess the assignment for 'suspect' contigs. It will rerun these steps a number of times in order to reach a convergence for reassigning haplotigs.

```
#!text
zz2_autoassign_contigs.pl  -s  coverage_stats.csv  -g  genome.fasta

REQUIRED:
-s      The coverage stats .csv file output from STEP 2
-g      The input assembly genome .fasta file (also needs the samtools index genome.fasta.fai file)

OPTIONAL:
-o      Output file name for the reassignment .tsv file. 
        DEFAULT = "suspect_contig_reassign.tsv"
-t      Threads to use for the blastn search and mummer alignments.
        DEFAULT = 4
-p      Max number of passes to perform, DEFAULT = 3. More than one purging
        pass is usually needed due to multiple overlapping haplotigs and repeat contigs.
-c      Prefix for the curated assembly, DEFAULT = "curated"

-m      Maxmatch cutoff percentage. Used to determine if a contig is a 
        repetitive sequence. DEFAULT = 250
-a      Bestmatch cutoff percentage. Used to determine if a contig is a
        haplotig. DEFAULT = 75

-u      Produce dotplots for unknown contigs only. DEFAULT is to produce
        dotplots for both assigned and unassigned.
```

### ALL DONE! 

You will have three files:

- <prefix>.fasta: These are the curated primary contigs
- <prefix>.haplotigs.fasta: These are all the haplotigs identified in the initial input assembly.
- <prefix>.artefacts.fasta: These are the low/high coverage contigs identified in STEP 2; most likely repeats and other assembly artefacts.

You can do some further reassigning by hand if you wish, see the steps below.


#### OPTIONAL MANUAL STEP - FURTHER CURATION

Go through `suspect_contig_reassign.tsv`, look at the corresponding dotplots for each unknown contig, and make your own assessment. Modify `suspect_contig_reassign.tsv` by hand by adding or changing reassignment keys.

```
#!text

Reassign_keys: h = haplotig, r = repeat/assembly junk, c = crop, ? = ¯\_(ツ)_/¯
Only 'h', 'r', and 'c' flagged contigs will be reassigned.
For cropping: positions are 1-indexed, can use "start" and "end" for start/end of contig
NOTE: The crop region you specify here is the region you wish to KEEP as a primary contig
(the rest of the contig will be output with the haplotigs).
SUSPECT_CONTIG    TOP_MATCH  SECOND_MATCH  MAXMATCHCOV      BESTMATCHCOV      REASSGIN_KEY  CROP_START  CROP_END
contig1            hit1       hit2          10.00            5.00              
contig2            hit1       hit2          95.00            93.00             h
contig3            hit1       hit2          400.00           98.00             r
contig4            hit1       hit2          50.00            50.00             c             50000       end
```

[Example: x-axis contig is a haplotig](https://bitbucket.org/mroachawri/purge_haplotigs/src/cf363f94c00fd865891a0469675d6df4a0813820/examples/example_haplotig.png)

[Example: x-axis contig is partially a haplotig](https://bitbucket.org/mroachawri/purge_haplotigs/src/cf363f94c00fd865891a0469675d6df4a0813820/examples/example_partial_haplotig.png)

[Example: x-axis contig is collapsed repeat/assembly junk](https://bitbucket.org/mroachawri/purge_haplotigs/src/cf363f94c00fd865891a0469675d6df4a0813820/examples/example_repetitive_junk_contig.png)

#### OPTIONAL - STEP 4

Run the fourth script to re-create your manually-reassigned curated assembly.

```
#!text

zz3_reassign_contigs.pl  -t suspect_contig_reassign.tsv  -g genome.fasta  -o output_prefix

-t/table            The output .tsv file from previous step (either manually edited after
                    reviewing the dotplot files, or unedited and using the auto-
                    assignments).

-g/genome           The curated assembly from STEP 3; NOT the original input assembly.

-o/output_prefix    Prefix for the output files. The following files will be created: 
                    <prefix>.fasta              - new primary contig assembly
                    <prefix>.haplotigs.fasta    - primary contigs reassigned as haplotigs
                    <prefix>.artefacts.fasta    - repeats and junk, probably not useful
```
#### DONE