#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PipeUtils;



#---SET UP LOGGING---
our $LOG;
my $TMP_DIR = "tmp_purge_haplotigs";
(-d $TMP_DIR) or mkdir $TMP_DIR;
open $LOG, ">", "$TMP_DIR/purge_haplotigs_readhist.log" or die "failed to open log file for writing";



# pre flight
if (check_programs("bedtools", "Rscript", "samtools")){
    msg("ALL DEPENDENCIES OK");
} else {
    err("ONE OR MORE DEPENDENCIES MISSING");
}
my $SCRIPT = "$Bin/../scripts/";



my $usage = "
USAGE:
purge_haplotigs  readhist  aligned.bam

REQUIRED:
aligned.bam   Samtools-indexed bam file of aligned and sorted reads/subreads to the reference

";



# parse and check file
my $bamfile = shift or die $usage;

if (!(check_files("$bamfile"))){
    die $usage;
}

# run bedtools genomecov if needed
if (!(-s "$TMP_DIR/$bamfile.genecov")){
    runcmd({ command => "bedtools genomecov -ibam $bamfile -max 200 > $TMP_DIR/$bamfile.gencov.temp 2> $TMP_DIR/bedtools.genomecov.stderr  &&  mv $TMP_DIR/$bamfile.gencov.temp $bamfile.gencov", logfile => "$TMP_DIR/bedtools.genomecov.stderr" });
} else {
    msg("$TMP_DIR/$bamfile.genecov found, skipping bedtools genomecov step");
}

# make the histogram csv
if (!(-s "$TMP_DIR/$bamfile.histogram.csv")){
    runcmd( { command => "grep genome $bamfile.gencov | awk '{ print \$2 \",\" \$3 }' > $TMP_DIR/$bamfile.histogram.csv" } );
} else {
    msg("$bamfile.histogram.csv found, skipping step");
}

# make the histogram png
if (!(-s "$bamfile.histogram.png")){
    runcmd( { command => "$SCRIPT/gen_histogram.Rscript $TMP_DIR/$bamfile.histogram.csv $bamfile.histogram.png 2> $TMP_DIR/gen_histogram.stderr", logfile => "$TMP_DIR/gen_histogram.stderr" });
} else {
    msg("$bamfile.histogram.png found, skipping Rscript step");
}

# done
msg("
purge_haplotigs readhist has finished! 

Check '$bamfile.histogram.png' to observe where your haploid and diploid peaks are
and choose your low, midpoint, and high cutoffs (check the example histogram png 
files in this git to give you an idea). You will need '$bamfile.gencov' and the 
cutoffs for the next step 'purge_haplotigs contigcov'
");

exit(0);

