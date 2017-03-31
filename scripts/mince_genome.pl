#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "
Usage:

get_seqs.pl  -file  <in.fasta>  -out  <DIR_NAME>

INPUT:
-file       Input fasta file

OUTPUT:
-out        Name of directory for output

";

my $in;
my $outDIR;
my $IN;
my $OUT;

GetOptions(
    "file=s" => \$in,
    "out=s" => \$outDIR
) or die $usage;

if ( !($in) || !($outDIR) ){
    die $usage;
}

open($IN, $in) or die "couldn't open $in for reading\n";

if ( !(-d $outDIR) ){
    mkdir $outDIR;
}

while (<$IN>) {
    if ($_ !~ /^>/) {
        print $OUT $_;
    } else {
        my $id;
        if ($_ =~ /^>([a-zA-Z0-9_-]+)[\s\|]/){
            $id = $1;
            close $OUT if ($OUT);
            open $OUT, ">", "$outDIR/$id.fasta" or die "failed to open \"$outDIR/$id.fasta\" for writing\n";
            print $OUT $_;
        } else {
            print STDERR "failed to get id from fasta file\n";
            exit(1);
        }
    }
}

close($OUT);
close($IN);
exit(0);

