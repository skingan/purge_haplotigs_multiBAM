#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PipeUtils;

my $usage = "
USAGE:
purge_haplotigs  readhist  <bam>

REQUIRED:
<bam>   Sorted bam file of reads aligned to your genome assembly

";

if (check_programs("bedtools", "Rscript")){
    msg("ALL DEPENDENCIES OK");
} else {
    err("ONE OR MORE DEPENDENCIES MISSING");
}

my $bamfile = shift or die $usage;
my $TEMP = "tmp_purge_haplotigs";
my $SCRIPT = "$Bin/../scripts/";

if (!(check_files("$bamfile"))){
    die $usage;
}

if (!(-e "$bamfile.gencov.done") || !(-s "$bamfile.genecov")){
    runcmd("bedtools genomecov -ibam $bamfile -max 200 > $bamfile.gencov");
    qruncmd("touch $bamfile.genecov.done");
} else {
    msg("$bamfile.genecov found, skipping bedtools genomecov step");
}

if (!(-d $TEMP)){
    mkdir $TEMP or err("failed to create temp folder $TEMP");
}

if (!(-s "$TEMP/$bamfile.histogram.csv")){
    runcmd("grep genome $bamfile.gencov | awk '{ print \$2 \",\" \$3 }' > $TEMP/$bamfile.histogram.csv");
} else {
    msg("$bamfile.histogram.csv found, skipping bash step");
}

if (!(-s "$bamfile.histogram.png")){
    runcmd("$SCRIPT/gen_histogram.Rscript $TEMP/$bamfile.histogram.csv $bamfile.histogram.png");
} else {
    msg("$bamfile.histogram.png found, skipping Rscript step");
}

msg("
purge_haplotigs readhist has finished! 

Check '$bamfile.histogram.png' to observe where your haploid and diploid peaks are
and choose your low, midpoint, and high cutoffs (check the example histogram png 
files in this git to give you an idea). You will need '$bamfile.gencov' and the 
cutoffs for the next step 'purge_haplotigs contigcov'
");
exit(0);

