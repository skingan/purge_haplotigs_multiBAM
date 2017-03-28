#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "
Usage:

returnseq.pl  -file  <in.fasta>  -list  <contigs.list>  -out  <DIR_NAME>

INPUT:
-file       Input fasta file

-list       List file of complete seq names to match

OUTPUT:
-out        Name of directory for output

";

my $in;
my $list;
my $outDIR;
my %contigs;
my $IN;
my $LIST;
my $OUT;

GetOptions(
    "file=s" => \$in,
    "list=s" => \$list,
    "out=s" => \$outDIR
) or die $usage;

if ( !($in) || !($list) || !($outDIR) ){
    die $usage;
}

open($IN, $in) or die "couldn't open $in for reading\n";

open($LIST, $list) or die "couldn't open $list for reading\n";
while(<$LIST>){
    $_ =~ s/\s//g;
    $contigs{$_} = 1;
}

if ( !(-d $outDIR) ){
    mkdir $outDIR;
}

my $m = 0;
while (<$IN>) {
    if ($_ !~ /^>/) {
        print $OUT $_ if ($m == 1);
    } else {
        my $id;
        if ($_ =~ /^>([a-zA-Z0-9_-]+)\s/){
            $id = $1;
        } else {
            die "failed to get id from fasta file\n";
        }
        if ($contigs{$id}){
            $m=1;
            close $OUT if ($OUT);
            open $OUT, ">", "$outDIR/$id.fasta" or die "failed to open \"$outDIR/$id.fasta\" for writing\n";
            print $OUT $_;
        } else {
            $m = 0;
        }
    }
}

