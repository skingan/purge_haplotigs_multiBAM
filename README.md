# purge_haplotigs

series of scripts to guide identifying haplotigs in a heterozygous diploid assembly (FALCON or FALCON-unzip)


### Dependencies

- bash
- bedtools
- Rscript with ggplot2
- perl 
- blast (blastn and makeblastdb)
- make (for mummer)


### Install

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
- install mdp_mummer (I modified mummer with some fixes and to make nicer looking dotplots by default, otherwise it's essentially just the mummer package)

```
cd /path/to/purge_haplotigs/dotplot_maker/MDP_MUMer
./install.sh
```

 - add the purge_haplotigs bin to your system PATH (**not** the dotplot_maker/bin--the scripts will find the relative path to this directory).

```
cd /path/to/purge_haplotigs/bin
$PATH=$PATH:$PWD
```


### Usage

map pacbio subreads or some decent short reads to your genome (we want a library that produces nice even coverage; we also want a randombest alignment for multimappers), and sort the bam. 

- generate a coverage histogram (either manually or run the script)

```
zz0_coverage_hsitogram.sh aligned.sorted.bam
```

eyeball the histogram, you should have two peaks, one for haploid level of coverage, the other for diploid level of coverage. choose cutoffs for low coverage, low point between the peaks, and high coverage (see example .png file)

- run script to analyse the coverage on a contig by contig basis, produces a contig stats .csv file with flagged contigs for further analysis

```
zz1_analyse_gencov.pl  -i genecov.out  -o stats.csv  -l 30  -m 80  -h 145  [ -j 80  -s 80 ]

-i      The output of bedtools genomecov
-o      Output file name (csv format)
-l      The read depth low cutoff (use the histogram to eyeball these cutoffs)
-h      The read depth high cutoff
-m      The low point between the haploid and diploid peaks

-j      Flag contig as \"j\" (junk) if this percentage or greater of the contig is 
            low/high coverage (default=80)
-s      Flag contig as \"s\" (suspected haplotig) if this percentage or more of the
            contig is haploid level of coverage (default=80)

```

- Run script to guess the assignment of the flagged contigs (produces `suspect_contig_reassign.tsv` )

```
zz2_assign_contigs.sh  stats.csv  genome.fasta

NOTE: make sure the genome is indexed with samtools faidx (need genome.fasta.fai)
contig names should only contain alphanumerics, underscores, or hyphens (you'll need to edit the regex in the script otherwise)
```

OPTIONAL: go through `suspect_contig_reassign.tsv` and the produced dotplots and make your own assessment; modify the .tsv file by hand. Otherwise the next step will just reassign the automatically identified contigs.

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

- run script to create your new, cleaner assembly

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

You can concatonate the reassigned <prefix>.haplotigs.fasta to your original haplotigs
file after this step (for instance if you've used FALCON-unzip).