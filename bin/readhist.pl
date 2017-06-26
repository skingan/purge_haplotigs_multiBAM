#!/usr/bin/env perl

use strict;
use warnings;
use Time::Piece;
use FindBin qw($Bin);

my $usage = "
USAGE:
purge_haplotigs  readhist  <bam>

REQUIRED:
<bam>   Sorted bam file of reads aligned to your genome assembly

";

if (chkprog("bedtools", "Rscript")){
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

if (!(-s "$bamfile.gencov")){
    runcmd("bedtools genomecov -ibam $bamfile -max 200 > $bamfile.gencov");
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
'purge_haplotigs readhist' has finished. Check '$bamfile.histogram.png' to observe where your haploid and diploid peaks are and choose your low, midpoint, and high cutoffs (check the example histogram png files in this git to give you an idea). You will need '$bamfile.gencov' and the cutoffs for the next step 'purge_haplotigs contigcov'\n");
exit(0);

#--- 

sub print_message {
    my $t = localtime;
    my $line = $t->dmy . " " . $t->hms . " @_\n";
    print STDERR $line;
}

sub msg {
    print_message("INFO: @_");
}

sub err {
    print_message("ERROR: @_\n\nPurge_Haplotigs has failed.\n");
    exit(1);
}

sub runcmd {
    print_message("RUNNING: @_");
    system("@_") == 0 or err("Failed to run @_\nCheck $_[-1] for possible causes.");
    print_message("FINISHED: @_")
}

sub check_files {
    my $check=1;
    foreach(@_){
        if (!(-s $_)){
            print_message("ERROR: file \"$_\" does not exist or is empty");
            $check=0;
        }
    }
    return $check;
}

sub chkprog {
    my $chk=1;
    foreach my $prog (@_){
        print_message("CHECKING $prog");
        my $notexists = `type $prog 2>&1 1>/dev/null || echo 1`;
        if ($notexists){
            print_message("ERROR: missing program $prog");
            $chk = 0;
        } else {
            msg("$prog OK");
        }
    }
    return $chk;
}
