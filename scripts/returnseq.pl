#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;


my $usage = "
Usage:

returnseq.pl  -f  <in.fasta>  ( -id 000000F  ||  -list contigs.fofn )  >  <out.fasta>

NEEDED:
-file       Input fasta file

NEED ONE OF:
-id         Complete seq name to match

-list       List file of complete seq names to match

";

my $in;
my $list;
my $ms;
my %contigs;
my $IN;
my $LIST;

GetOptions(
    "file=s" => \$in,
    "id=s" => \$ms,
    "list=s" => \$list
) or die $usage;

if ( (!($ms) && !($list)) || ($ms) && ($list) ){
    die $usage;
}

open($IN, $in) or die "couldn't open $in for reading\n";

if ($list){
    print STDERR "Matching contigs:\n";
    open($LIST, $list) or die "couldn't open $list for reading\n";
    while(<$LIST>){
        $_ =~ s/\s//g;
        $contigs{$_} = 1;
        print STDERR "$_\n";
    }
}

my $m = 0;
while (<$IN>) {
    if ($_ !~ /^>/) {
        print STDOUT $_ if ($m == 1);
    } else {
        if ($list){
            my $id;
            if ($_ =~ /^>([a-zA-Z0-9_-]+)\s/){
                $id = $1;
            } else {
                die "failed to get id from fasta file\n";
            }
            if ($contigs{$id}){
                $m=1;
                print STDOUT $_;
            } else {
                $m = 0;
            }
        } else {
            if ($_ =~ /$ms\s/) {
                $m = 1;
                print STDOUT $_;
            } else {
                $m = 0;
            }
        }
    }
}

