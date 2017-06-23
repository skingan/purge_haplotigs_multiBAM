# purge_haplotigs

Pipeline to help in the curation of a heterozygous diploid genome assembly (such as using FALCON or FALCON-unzip).

#### The problem

Some parts of a genome may have a very high degree of heterozygosity. This causes contigs for both haplotypes of that part of the genome to be assembled as separate primary contigs, rather than as a contig and an associated haplotig. This can be an issue for downstream analysis whether you're working on the haploid or phased-diploid assembly.

#### The solution

Identify primary contigs that are haplotigs of other primary contigs, and add them to the haplotigs. The pipeline will use mapped read coverage and blast/lastz alignments to determine which contigs to keep as primary contigs in order to create a highly contiguous haploid assembly. Dotplots are produced for all flagged contig matches, juxtaposed with read-coverage, to help the user determine the proper assignment of any remaining ambiguous contigs.

## Dependencies

- Bash
- BEDTools and SAMTools
- blast (blastn and makeblastdb)
- lastz - [website](http://www.bx.psu.edu/~rsharris/lastz/), [github](https://github.com/lastz/lastz)
- Perl 
- Rscript (with ggplot2 and scales)
- GNU Parallel (optional but highly recommended)


## Install (Ubuntu/Debian)

- Install dependencies, make sure they're in the system PATH
- Pull/clone this git
- Add the purge_haplotigs bin to your system PATH

```
#!bash
cd /path/to/purge_haplotigs/bin
PATH=$PATH:$PWD

# optionally add a line to your .bashrc
printf "\nPATH=\$PATH:$PWD\n" >> $HOME/.bashrc
```

That's it!

## RUNNING THE PIPELINE

#### PREPARATION

Map your PacBio subreads, or some decent long reads (or even short reads) to your genome assembly (use a library that produces an even coverage and use a 'randombest' alignment for multimappers). NOTE: you can use the phased-diploid assembly (primary contigs + haplotigs) or just the haploid assembly (primary contigs only) for this step. Sort and index the bam with SamTools. Index your genome.fasta file with samtools faidx if you haven't already done so.

#### STEP 1

Generate a coverage histogram by running the first script. You will also need the bedtools genecov output file for STEP 2.

```
#!text
zz0_coverage_hsitogram.sh aligned.bam
```

#### MANUAL STEP

You should have a bimodal histogram--one peak for haploid level of coverage, one peak for diploid level of coverage. NOTE: If you're using the phased assembly the diploid peak may be very small. Choose cutoffs for low coverage, low point between the peaks, and high coverage.

**Example histograms for choosing cutoffs:**

[PacBio subreads on Diploid-phased assembly (Primary + Haplotigs)](https://bitbucket.org/mroachawri/purge_haplotigs/src/16c0b3f8d1c7eaf32d2cac5441c6864a1cc92bd3/examples/phased_coverage_histogram.png)

[Illumina PE reads on Haploid assembly (Primary contigs)](https://bitbucket.org/mroachawri/purge_haplotigs/src/16c0b3f8d1c7eaf32d2cac5441c6864a1cc92bd3/examples/coverage_histogram.png)

#### STEP 2

Run the second script using the cutoffs from the previous step to analyse the coverage on a contig by contig basis. This script produces a contig `coverage_stats.csv` file with suspect contigs flagged for further analysis or removal depending on read-coverage.

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
-j      Auto-assign contig as "j" (junk) if this percentage or greater of the contig is
            low/high coverage (default=80, >100 = off)
-s      Auto-assign contig as "s" (suspected haplotig) if this percentage or less of the
            contig is diploid level of coverage (default=80)

```

#### STEP 3

Run the purging pipeline. This script will automatically run a BEDTools windowed coverage analysis, a blast search, and perform lastz alignments etc. to assess which contigs to reassign and which to keep. The pipeline will make several iterations of purging. If you rerun this step it won't by default rerun the BEDTools or Blastn analysis.

```
#!text
purge_haplotigs.pl  -genome genome.fasta  -coverage coverage_stats.csv -bam aligned.sorted.bam

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

### ALL DONE! 

#### You will have five files:

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

#### You will also get two directories of dotplots juxtaposed with read-coverage tracks:

- **dotplots_unassigned_contigs:** These are the dotplots for the suspect contigs that remain unassigned.
- **dotplots_reassigned_contigs:** These are the dotplots for the contigs that were reassigned. 

### FURTHER CURATION

You can go through the dotplots and check the assignments. You might wish to move reassigned contigs back into the primary contig pool or purge contigs that were too ambiguous for automatic reassignment. Below are some examples of what might occur. 

[A haplotig](https://bitbucket.org/mroachawri/purge_haplotigs/src/16c0b3f8d1c7eaf32d2cac5441c6864a1cc92bd3/examples/haplotig.png)

[Contig is mostly haplotig](https://bitbucket.org/mroachawri/purge_haplotigs/src/16c0b3f8d1c7eaf32d2cac5441c6864a1cc92bd3/examples/haploid_diploid_hemizygous.png) - This example has part of the contig with a diploid level of coverage and part at haploid level with poor alignment to either reference (possibly hemizygous region).

[Haplotig with many tandem repeats](https://bitbucket.org/mroachawri/purge_haplotigs/src/16c0b3f8d1c7eaf32d2cac5441c6864a1cc92bd3/examples/repeat_rich.png)

[Haplotig is palindrome](https://bitbucket.org/mroachawri/purge_haplotigs/src/16c0b3f8d1c7eaf32d2cac5441c6864a1cc92bd3/examples/haplotig_with_palindrome.png)

[Contig circling through string graph 'knots'](https://bitbucket.org/mroachawri/purge_haplotigs/src/16c0b3f8d1c7eaf32d2cac5441c6864a1cc92bd3/examples/repeats_string_graph_short_cut.png) - while this may be a valid string graph path it will likely still confound short read mapping to the haploid assembly.
