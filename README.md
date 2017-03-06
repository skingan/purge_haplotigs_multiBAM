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

Generate a coverage histogram (either by running the first script or using your favourite method)

```
#!text
zz0_coverage_hsitogram.sh aligned.sorted.bam
```

#### MANUAL STEP

eyeball the histogram, you should have two peaks, one for haploid level of coverage, the other for diploid level of coverage. Choose cutoffs for low coverage, low point between the peaks, and high coverage.

[Example histogram and choosing cutoffs](https://bitbucket.org/mroachawri/purge_haplotigs/src/cf363f94c00fd865891a0469675d6df4a0813820/examples/example_histogram.png)

#### STEP 2

Run the second script to analyse the coverage on a contig by contig basis--produces a contig `stats.csv` file with flagged contigs for further analysis

```
#!text
zz1_analyse_gencov.pl  -i genecov.out  -o stats.csv  -l 30  -m 80  -h 145  [ -j 80  -s 80 ]

-i      The output of bedtools genomecov
-o      Output file name (csv format)
-l      The read depth low cutoff (use the histogram to eyeball these cutoffs)
-h      The read depth high cutoff
-m      The low point between the haploid and diploid peaks

-j      Flag contig as "j" (junk) if this percentage or greater of the contig is 
            low/high coverage (default=80)
-s      Flag contig as "s" (suspected haplotig) if this percentage or more of the
            contig is haploid level of coverage (default=80)

```

#### STEP 3

Run the third script to guess the assignment of the flagged contigs (produces `suspect_contig_reassign.tsv` )

```
#!text
zz2_assign_contigs.sh  stats.csv  genome.fasta
```

#### OPTIONAL MANUAL STEP

Go through `suspect_contig_reassign.tsv`, look at the corresponding dotplots for each unknown contig, and make your own assessment. Modify `suspect_contig_reassign.tsv` by hand. Otherwise the next step will reassign the automatically identified contigs only.

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

#### STEP 4

Run the fourth script to create your new, cleaner assembly

```
#!text

zz3_reassign_contigs.pl  -t suspect_contig_reassign.tsv  -g genome.fasta  -o output_prefix

-t/table            The output .tsv file from previous step (either manually edited after
                    reviewing the dotplot files, or unedited and using the auto-
                    assignments).

-g/genome           The original genome .fasta (primary contigs from FALCON/FALCON-unzip).

-o/output_prefix    Prefix for the output files. The following files will be created: 
                    <prefix>.fasta              - new primary contig assembly
                    <prefix>.haplotigs.fasta    - primary contigs reassigned as haplotigs
                    <prefix>.artefacts.fasta    - repeats and junk, probably not useful

-force              Force reassignment (if contig is flagged for reassigning but is also
                    a reference for flagging another contig) off=script will try to 
                    determine which to keep.


```
#### DONE

You can concatenate the reassigned <prefix>.haplotigs.fasta to your original haplotigs
file after this step (for instance if you've used FALCON-unzip).