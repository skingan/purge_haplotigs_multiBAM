#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
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
my $SCRIPT = "$RealBin/../scripts/";



my $usage = "
USAGE:
purge_haplotigs  readhist  aligned.bam

REQUIRED:
aligned.bam   BAM file of aligned and sorted reads/subreads to the reference

";



# parse and check file
my $bamfile = shift or die $usage;
my $bamfilename = $bamfile;
$bamfilename =~ s/.*\///;

if (!(check_files("$bamfile"))){
    die $usage;
}


# FILENAMES FOR CLARITY
my $gencov_file = "$bamfilename.genecov";
my $csv_file = "$TMP_DIR/$bamfilename.histogram.csv";
my $png_file = "$bamfilename.histogram.png";



# run bedtools genomecov if needed
if (!(-s $gencov_file)){
    runcmd({ command => "bedtools genomecov -ibam $bamfile -max 200 > $TMP_DIR/$gencov_file.temp 2> $TMP_DIR/bedtools.genomecov.stderr  &&  mv $TMP_DIR/$gencov_file.temp $gencov_file", logfile => "$TMP_DIR/bedtools.genomecov.stderr" });
} else {
    msg("$gencov_file found, skipping bedtools genomecov step");
}

# make the histogram csv
if (!(-s $csv_file)){
    runcmd( { command => "grep genome $gencov_file | awk '{ print \$2 \",\" \$3 }' > $csv_file" } );
} else {
    msg("$csv_file found, skipping step");
}

# make the histogram png
if (!(-s $png_file)){
    runcmd( { command => "$SCRIPT/gen_histogram.Rscript $csv_file $png_file 2> $TMP_DIR/gen_histogram.stderr", logfile => "$TMP_DIR/gen_histogram.stderr" });
} else {
    msg("$png_file found, skipping Rscript step");
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

