#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);

my $usage = "
Usage:
analyze_blastn_output.pl  -f genome.fasta  -b blastn_out.fmt6.gz  -t CONTIG_REASSIGN.tsv  [ -m 150  -a 90 ]

-f      genome.fasta file, needs to be indexed with samtools faidx
-b      output from previous pipeline step, blastn, -outfmt 6, optionally gzipped
-t      output .tsv file which then needs to be edited and fed into next step

OPTIONAL:
-m      maxmatch cutoff - for finding collapsed repeat contigs. -m 150 means that if a contig 
            matches 150 % or more to another contig then it's a collased repeat/assembly junk.
-a      best alignment cutoff - for identifying haplotigs. -a 90 means that if a contig is not
            a collapsed repeat, and maps 90 % or more to another contig, then it's a haplotig.

";

# args
my $fasta;
my $blast;
my $tblout;
my $max_match_cutoff = 150;
my $align_match_cutoff = 75;

GetOptions(
    "fastafai=s" => \$fasta,
    "blast=s" => \$blast,
    "table=s" => \$tblout, 
    "maxmatch=s" => \$max_match_cutoff,
    "alignmatch" => \$align_match_cutoff
) or die $usage;

if (!($fasta) || !($blast) || !($tblout) ){
    die $usage;
}

# genome index
my $fai = "$fasta.fai";
if (!(-e $fai)){
    die "$fasta needs to be indexed with samtools faidx; cant find $fasta.fai\n";
}

# tempfolder
my $temp = "tmp_blastn_analysis";
if (!(-d $temp)){
    mkdir $temp;
}

# dotplot folder
my $dotunk = "dotplots_unknowns";
my $dotcall = "dotplots_assigned";
if (!(-d $dotunk)){
    mkdir $dotunk;
}
if (!(-d $dotcall)){
    mkdir $dotcall;
}

# filehandles
my $FAI;
my $BLST;
my $OUT;

# global vars
my %length; # $length{contig} = 2000000
my %hits; # $hits{contig}{1} = contig; $hits{contig}{2} = anothercontig
my %RBHs; # $RBHs{largercontig} = smallercontig
my %purged; # $purged{contig} = 1

my $MDPbin = "$Bin/../dotplot_maker/bin/";

# open files
open($FAI, $fai) or die "failed to open $fai for reading\n";
if ($blast =~ /\.gz$/){
    open($BLST, "gunzip -c $blast |") or die "failed to open gunzip pipe from $blast\n";
} else {
    open($BLST, $blast) or die "failed to open file $blast for reading\n";
}
open($OUT, ">$tblout") or die "failed to open $tblout for writing\n";

# read in the contig lengths
while (<$FAI>){
    my @line = split(/\s+/, $_);
    $length{$line[0]} = $line[1];
}
close($FAI);

# read in blast output file
while (<$BLST>){
    my @line = split(/\s+/, $_);
    if (!($hits{$line[0]}{1})){
        $hits{$line[0]}{1} = $line[1] if ($line[0] ne $line[1]); # shouldn't be self hits but just in case
    } elsif (!($hits{$line[0]}{2}) && ($hits{$line[0]}{1} ne $line[1])){
        $hits{$line[0]}{2} = $line[1] if ($line[0] ne $line[1]);
    }
}
close($BLST);

# characterise hits
foreach my $contig (keys(%hits)){
    # find reciprocal best hits, this will only happen when both contigs were flagged as suspected haplotigs, which
        # hopefully should be the case for most of the problematic contigs
    if ( ($hits{$contig}{1}) && ($hits{$hits{$contig}{1}}{1}) ){
        if ( $hits{$hits{$contig}{1}}{1} eq $contig ){
            if ( $length{$contig} > $length{$hits{$contig}{1}} ){
                $RBHs{$contig} = $hits{$contig}{1};
            } else {
                $RBHs{$hits{$contig}{1}} = $contig;
            }
        }
    }
    
    # ignore if contig is larger that both its hit seqs
    elsif ( $length{$contig} > ($length{$hits{$contig}{1}} + $length{$hits{$contig}{2}}) ){
        $purged{$contig} = 1;
    }
}

# TESTING
print $OUT "#SUSPECT_CONTIG\tTOP_MATCH\tSECOND_MATCH\tMaxMatchCov\tBestMatchCov\tREASSGIN_KEY\tCROP_START\tCROP_END\n";
print $OUT "#Reciprocal_best_hits\n";
foreach my $contig (keys(%RBHs)){
    if (!($purged{$RBHs{$contig}})){
        # we assume the smaller contig is the haplotig
        print STDERR "$Bin/returnseq.pl -f $fasta -i $hits{$RBHs{$contig}}{1} > $temp/ref.fasta\n";
        `$Bin/returnseq.pl -f $fasta -i $hits{$RBHs{$contig}}{1} > $temp/ref.fasta`;
        print STDERR "$Bin/returnseq.pl -f $fasta -i $hits{$RBHs{$contig}}{2} >> $temp/ref.fasta\n";
        `$Bin/returnseq.pl -f $fasta -i $hits{$RBHs{$contig}}{2} >> $temp/ref.fasta`;
        my $assign = guess_assignment($RBHs{$contig});
        print $OUT "$RBHs{$contig}\t$contig\t$hits{$RBHs{$contig}}{2}\t$assign\n";
    }
}
print $OUT "#Everything_else\n";
foreach my $contig (sort(keys(%hits))){
    if (!($purged{$contig})){
        if ($hits{$contig}{1}){
            `$Bin/returnseq.pl -f $fasta -i $hits{$contig}{1} > $temp/ref.fasta`;
            `$Bin/returnseq.pl -f $fasta -i $hits{$contig}{2} >> $temp/ref.fasta`;
            my $assign = guess_assignment($contig);
            print $OUT "$contig\t$hits{$contig}{1}\t$hits{$contig}{2}\t$assign\n";
        }
    }
}
close($OUT);

print STDERR "\nTable is saved in the following format (check and update the table, feed into next step of pipeline)\n";
print STDERR "Reassign_keys: h = haplotig, r = repeat/assembly junk, c = crop, blank = no reassignment, ? = ¯\\_(^_^)_/¯\n";
print STDERR "For cropping: positions are 1-indexed, can use \"start\" and \"end\" for start/end of contig\n\n";
print STDERR "The crop region you specify here is the region you wish to KEEP.\n\n";
print STDERR "#SUSPECT_CONTIG    TOP_MATCH  SECOND_MATCH  MAXMATCHCOV      BESTMATCHCOV      REASSGIN_KEY  CROP_START  CROP_END\n";
print STDERR "#contig1           hit1       hit2          10.00            5.00              \n";
print STDERR "#contig2           hit1       hit2          95.00            93.00             h\n";
print STDERR "#contig3           hit1       hit2          400.00           98.00             r\n";
print STDERR "#contig4           hit1       hit2          50.00            50.00             c             50000       end\n\n";

exit(0);

sub guess_assignment{
    my $query = $_[0];
    my $ref = "$temp/ref.fasta";
    # get seq
    print STDERR "$Bin/returnseq.pl -f $fasta -i $query > $temp/tmp_query.fasta\n";
    `$Bin/returnseq.pl -f $fasta -i $query > $temp/tmp_query.fasta`;
    
    # nucmer
    print STDERR "$MDPbin/MDP_nucmer --maxmatch -p $temp/tmp $temp/tmp_query.fasta $ref\n";
    `$MDPbin/MDP_nucmer --maxmatch -p $temp/tmp $temp/tmp_query.fasta $ref`;
    
    # deltas
    print STDERR "$MDPbin/MDP_delta-filter -m $temp/tmp.delta > $temp/tmp.m.delta\n";
    `$MDPbin/MDP_delta-filter -m $temp/tmp.delta > $temp/tmp.m.delta`;
    print STDERR "$MDPbin/MDP_delta-filter -1 $temp/tmp.delta > $temp/tmp.1.delta\n";
    `$MDPbin/MDP_delta-filter -1 $temp/tmp.delta > $temp/tmp.1.delta`;
    
    # all-match coverage
    print STDERR "$MDPbin/MDP_show-coords -b -c $temp/tmp.m.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'\n";
    my $maxmatch = `$MDPbin/MDP_show-coords -b -c $temp/tmp.m.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'`;
    $maxmatch =~ s/\s//g;
    print STDERR "MAXMATCH coverage = $maxmatch\n";
    
    # best-align coverage
    print STDERR "$MDPbin/MDP_show-coords -b -c $temp/tmp.1.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'\n";
    my $alignmatch = `$MDPbin/MDP_show-coords -b -c $temp/tmp.1.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'`;
    $alignmatch =~ s/\s//g;
    print STDERR "BESTMATCH coverage = $alignmatch\n";
    
    my $a = "$maxmatch\t$alignmatch\t";
    
    # guess the assignment and make dotplots
    if ( ($maxmatch >= $max_match_cutoff) && ($alignmatch >= $align_match_cutoff) ) {
        $purged{$query}=1;
        $a .= "r";
        print STDERR "\n### $query = repeat/assembly junk\n\n";
        print STDERR "$MDPbin/MDP_mummerplot --fat -p $dotcall/$query $temp/tmp.m.delta\n";
        `$MDPbin/MDP_mummerplot --fat -p $dotcall/$query $temp/tmp.m.delta`;
        unlink "$dotcall/$query.filter";
        unlink "$dotcall/$query.fplot";
        unlink "$dotcall/$query.rplot";
        unlink "$dotcall/$query.gp";
        
    } elsif ($alignmatch >= $align_match_cutoff){
        $purged{$query}=1;
        $a .= "h";
        print STDERR "\n### $query = haplotig\n\n";
        print STDERR "$MDPbin/MDP_mummerplot --fat -p $dotcall/$query $temp/tmp.m.delta\n";
        `$MDPbin/MDP_mummerplot --fat -p $dotcall/$query $temp/tmp.m.delta`;
        unlink "$dotcall/$query.filter";
        unlink "$dotcall/$query.fplot";
        unlink "$dotcall/$query.rplot";
        unlink "$dotcall/$query.gp";
        
    } else {
        $purged{$query}=1;
        $a .= "?";
        print STDERR "\n### $query = unsure, use your mk-I eyeballs\n\n";
        `$MDPbin/MDP_mummerplot --fat -p $dotunk/$query $temp/tmp.m.delta`;
        unlink "$dotunk/$query.filter";
        unlink "$dotunk/$query.fplot";
        unlink "$dotunk/$query.rplot";
        unlink "$dotunk/$query.gp";
    }
    
    # clean-up
    print STDERR "rm $temp/tmp_query.fasta $temp/tmp.delta $temp/tmp.m.delta $temp/tmp.1.delta $temp/ref.fasta\n";    
    unlink "$temp/tmp_query.fasta";
    unlink "$temp/tmp.delta";
    unlink "$temp/tmp.m.delta";
    unlink "$temp/tmp.1.delta";
    unlink "$temp/ref.fasta";
    
    return($a);
}


