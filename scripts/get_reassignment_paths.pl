#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;



#---PARAMS---

my $reassignments_tsv;
my $reassignment_path_tsv;

my $usage = "
Usage:
get_reassignment_paths.pl  -i curated.reassignments.tsv  -o curated.reassignment_paths.tsv

-i      input, the concatonated reassignments tsv file from zz2_autoassign_contigs.pl

-o      output, tsv file of the 'paths' showing the relationships between the contigs 
        and the order of removal.

";



#---OPTIONS---

GetOptions (
    "in=s" => \$reassignments_tsv,
    "out=s" => \$reassignment_path_tsv
) or die $usage;

if (!($reassignments_tsv) || !($reassignment_path_tsv)){
    die "missing input\n$usage";
}



#---FILEHANDLES---

my $IN;
my $OUT;



#---GLOBAL PARAMS---

my %list;       # $list{reassigned_contig}{r} = ref_contig
                #                         {a} = assignment

my %ref;        # @($ref{ref_contig}) = (reassigned_contig1, reassigned_contig2, ...)

my @primaries;          # contigs that were not reassigned
my @current_path;       # ("contig,PRIMARY -> ", "contig,HAPLOTIG -> ", "contig,REPEAT -> ", etc...
my $current_depth;      # depth of recursive subroutine, used to modify @current_path

#---PROGRAM---

# open files
open $IN, $reassignments_tsv or die "Failed to open $reassignments_tsv for reading";
open $OUT, ">", $reassignment_path_tsv or die "Failed to open $reassignment_path_tsv for writing";

# read in TSV
while(<$IN>){
    next if ($_ =~ /^#/);
    my @line = split(/\t/, $_);
    next if ($line[1] eq "NA");
    
    if ($list{$line[0]}){
        die "Error, contig $line[0] appears reassigned multiple time in TSV\n";
    }
    
    $list{$line[0]}{r} = $line[1];
    $list{$line[0]}{a} = $line[2];
    push @{$ref{$line[1]}}, $line[0];
}

# get list of contigs that are references and are not reassigned
foreach my $contig (sort(keys(%ref))){
    if (!($list{$contig})){
        push @primaries, $contig;
    }
}

# now to remember how to write a recursive subroutine...
foreach my $ctg (@primaries){
    @current_path = ("$ctg,PRIMARY");
    $current_depth = 0;
    find_path($ctg);
    print $OUT "\n";
}



sub find_path {
    my $ref_ctg = $_[0];
    
    if (!($ref{$ref_ctg})){
        print $OUT @current_path;
        print $OUT "\n";
    } else {
        foreach my $ctg (@{$ref{$ref_ctg}}){
            push @current_path, " <- $ctg,$list{$ctg}{a}";
            $current_depth += 1;
            
            find_path($ctg);
            
            pop @current_path;
            $current_depth -= 1;
            for (my $i=0; $i<=$current_depth; $i++){
                $current_path[$i] =~ s/./ /g;
            }
            $current_path[$current_depth] =~ s/ $/|/;
        }
    }
}









