#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "
Usage:
reassign_contigs.pl  -t contig_table.tsv  -g genome.fasta  -o output_prefix

-t/table            The output .tsv file from previous step (either manually edited after
                    reviewing the dotplot files, or unedited and using the auto-
                    assignments).
-g/genome           The original genome .fasta (primary contigs from FALCON/FALCON-unzip).
-o/output_prefix    Prefix for the output files. The following files will be created: 
                    <prefix>.fasta              - new primary contig assembly
                    <prefix>.haplotigs.fasta    - primary contigs reassigned as haplotigs
                    <prefix>.artefacts.fasta    - repeats and very low/high coverage, 
                                                  contigs, probably not useful
-force              Force reassignment (if contig is flagged for reassigning but is also
                    a reference for flagging another contig) off=script will try to 
                    determine which to keep.

You can concatonate the reassigned <prefix>.haplotigs.fasta to your original haplotigs
file after this step (for instance if you've used FALCON-unzip).

";

# args
my $table;
my $genome;
my $outprefix;
my $force;

GetOptions (
    "table=s" => \$table,
    "genome=s" => \$genome,
    "output_prefix=s" => \$outprefix,
    "force" => \$force
) or die $usage;

if ( !($table) || !($genome) || !($outprefix) ){
    die $usage;
}

# filehandles
my $OG;
my $OH;
my $OA;
my $IT;
my $IG;
my $JL;
my $FF;

# outfilenames
my $outG = "$outprefix.fasta";
my $outH = "$outprefix.haplotigs.fasta";
my $outA = "$outprefix.artefacts.fasta";

# global variables
my %table;      # $table{"contig"}{1} = "contig-id"        top 2 hits for contig
                #                 {2} = "contig-id"
                #                 {A} = "r/h/c"
                #                 {M} = 200.00             maxmatch coverage
                #                 {B} = 95.00              bestmatch coverage
                #                 {C1} = INT               crop start
                #                 {C2} = INT               crop_stop
my @contigs_for_reassigning;

# other
my $junklist = "tmp_reassign_contigs/junk.list";
my %junkcontigs;
my %contig_length;

if ( (-s $outG) || (-s $outH) || (-s $outA) ){
    die "One or more of  the outifles:\n$outG\n$outH\n$outA\nalready exists, exiting...\n";
}

# fasta .fai file
my $fastafai = "$genome.fai";
if (!(-s $fastafai)){
    die "genome index $fastafai appears to be missing, have you indexed with samtools faidx?\n";
}

# open files
open($IT, $table) or die "failed to open $table for reading\n";
open($IG, $genome) or die "failed to open $genome for reading\n";
open($JL, $junklist) or die "failed to open $junklist for reading\n";
open($FF, $fastafai) or die "failed to open $fastafai for reading\n";
open($OG, ">$outG") or die "failed to open $outG for writing\n";
open($OH, ">$outH") or die "failed to open $outH for writing\n";
open($OA, ">$outA") or die "failed to open $outA for writing\n";


# read in the list of junk contigs
($_ =~ s/\s//g, $junkcontigs{$_} = 1) while(<$JL>);
close($JL);

# read in contig lengths
while(<$FF>){
    my @line = split(/\s/, $_);
    $line[0] =~ s/\|.+$//;
    $contig_length{$line[0]} = $line[1];
}
close($FF);

# read in table file
while(<$IT>){
    # skip comments
    next if ($_ =~ /^#/);
    
    my @line = split(/\s/, $_);
    if (scalar(@line) < 6){
        err("ERROR: Columns missing in table .tsv file at line:\n@line\nexiting...");
    }
    
    # skip if nothing to be done with contig
    if ($line[5] !~ /^[rRhHcC]$/){
        print STDERR "INFO: nothing to be done for $line[0], skipping\n";
        next;
    }
    
    push @contigs_for_reassigning, $line[0];
    $table{$line[0]}{1} = $line[1];
    $table{$line[0]}{2} = $line[2];
    $table{$line[0]}{"A"} = $line[5];
    $table{$line[0]}{"A"} =~ tr/RHC/rhc/;
    $table{$line[0]}{"M"} = $line[3];
    $table{$line[0]}{"B"} = $line[4];
    if ($line[5] =~ /c/i){
        $table{$line[0]}{"C1"} = $line[6];
        $table{$line[0]}{"C2"} = $line[7];
        (print STDERR "$_ , ") foreach (@line);
        print STDERR "\n";
    }
}
close($IT);

# identify if any reassigned contigs are references for reassignment - looking at best hits only
foreach my $contig (@contigs_for_reassigning){
    if ($table{$table{$contig}{1}}){
        if ( ($table{$table{$contig}{1}}{"A"} =~ /[hrc]/) && ($table{$contig}{"A"} =~ /[hrc]/) ){
            print STDERR "INFO: $contig is both flagged for reassignment AND is a reference for flagging another contig for reassingmnet\n";
            pick_best_reference($contig) if (!($force));
        }
    }
}

# read through genome, reassign and print to files on-the-fly
my $current_contig;
my $current_seq;
ITR: while(<$IG>){
    if ($_ =~ /^>/){
        my $id;
        if ($_ =~ /^>([a-zA-Z0-9-_]+)[\s\|]/){  
            $id = $1;
        } else {
            err("ERROR: illegal characters in seq ID?\n$_\n");
        }

        if (!($current_contig)){
            $current_contig = $id;
            next ITR;
        }
        # print to junk fasta if junk contig, otherwise process seq
        if ($junkcontigs{$current_contig}){
            print_seq($OA, "$current_contig\_REPEAT/JUNK", $current_seq);
        } else {        
            write_current_seq();
        }
        
        # reset
        $current_contig = $id;
        $current_seq = "";
        
    } else {
        $_ =~ s/\s//g;
        $current_seq .= $_;
    }
}
# process the final seq
write_current_seq();

exit(0);



### SUBROUTINES ###

sub pick_best_reference{
    my $contig = $_[0];
    my $m_ctg = $table{$contig}{1};
    
    # if both flagged as haplotigs or crop
        # keep the longest contig
    if ( (($table{$contig}{"A"} eq "h") && ($table{$m_ctg}{"A"} eq "h")) ||
         (($table{$contig}{"A"} eq "c") && ($table{$m_ctg}{"A"} eq "c")) ){
        if ($contig_length{$contig} < $contig_length{$m_ctg}){
            $table{$m_ctg}{"A"} = "?";
            print STDERR "INFO:    keeping $m_ctg, reassigning $contig\n";
        } else {
            $table{$contig}{"A"} = "?";
            print STDERR "INFO:    keeping $contig, reassigning $m_ctg\n";
        }
    }
    
    # else if only one is repetitive
        # choose it
    elsif ( ( ($table{$contig}{"A"} =~ /[hc]/) && ($table{$m_ctg}{"A"} eq "r") ) || 
            ( ($table{$contig}{"A"} eq "r")    && ($table{$m_ctg}{"A"} =~ /[hc]/) ) ){
        if ($table{$contig}{"A"} =~ /[hc]/){
            $table{$contig}{"A"} = "?";
            print STDERR "INFO:    keeping $contig, reassigning $m_ctg\n";
        } else {
            $table{$m_ctg}{"A"} = "?";
            print STDERR "INFO:    keeping $m_ctg, reassigning $contig\n";
        }
    }
    
    # else if one crop one haplotig
        # choose haplotig
    elsif ( ( ($table{$contig}{"A"} eq "h") && ($table{$m_ctg}{"A"} eq "c") ) || 
            ( ($table{$contig}{"A"} eq "c") && ($table{$m_ctg}{"A"} eq "h") ) ){
        if ($table{$contig}{"A"} eq "c"){
            $table{$contig}{"A"} = "?";
            print STDERR "INFO:    keeping $contig, reassigning $m_ctg\n";
        } else {
            $table{$m_ctg}{"A"} = "?";
            print STDERR "INFO:    keeping $m_ctg, reassigning $contig\n";
        }
    }
    
    # else if BOTH repetitive or crop 
        # choose longest contig
    elsif ( ($table{$contig}{"A"} eq "r") && ($table{$m_ctg}{"A"} eq "r") ){
        if ($contig_length{$contig} < $contig_length{$m_ctg}){
            $table{$m_ctg}{"A"} = "?";
            print STDERR "INFO:    keeping $m_ctg, reassigning $contig\n";
        } else {
            $table{$contig}{"A"} = "?";
            print STDERR "INFO:    keeping $contig, reassigning $m_ctg\n";
        }
    }
    
    else {
        err("ERROR: unknown combination of reassigning flags\n$contig $table{$contig}{'A'}\n$m_ctg $table{$m_ctg}{'A'}\n");
    }
}

sub write_current_seq{
    if ($table{$current_contig}{"A"}){
        if ($table{$current_contig}{"A"} eq "r"){
            # print to junk file
            print_seq($OA, "$current_contig\_REPEAT/JUNK", $current_seq);
        } elsif ($table{$current_contig}{"A"} eq "h"){
            # print to haplotig file
            print_seq($OH, "$current_contig\_HAPLOTIG", $current_seq);
        } elsif ($table{$current_contig}{"A"} eq "c"){
            # crop and print to genome/haplotig
            crop_and_print();
        } else {
            # print to genome, nothing to be done with seq
            print_seq($OG, $current_contig, $current_seq);
        }
    } else {
        # print to genome, nothing to be done with seq
        print_seq($OG, $current_contig, $current_seq);
    }
}

sub crop_and_print{
    my $cropstart = $table{$current_contig}{"C1"};
    my $cropstop =  $table{$current_contig}{"C2"};
    
    if ($cropstart =~ /start/i){
        $cropstart = 0;
    }
    if ($cropstart <= 1){
        $cropstart = 0;
    }
    if ($cropstop =~ /end/i){
        $cropstop = (length($current_seq) - 1);
    }
    if ($cropstop >= length($current_seq)){
        $cropstop = (length($current_seq) - 1);
    }
    
    # integer checks
    if ($cropstart !~ /^\d+$/){
        err("ERROR: crop starting position not an integer\n$current_contig\n");
    }
    if ($cropstop !~ /^\d+$/){
        err("ERROR: crop stop position not an integer\n$current_contig\n");
    }
    
    # other checks
    if ($cropstart > $cropstop){
        err("ERROR: crop start position greater than crop stop position\n$current_contig\n");
    }
    
    # cut and print the seqs
        # middle crop
    if ( ($cropstart > 0) && ($cropstop < (length($current_seq) - 1)) ){
        my $pseq = substr($current_seq, $cropstart, ($cropstop - $cropstart));
        my $hseq1 = substr($current_seq, 0, $cropstart);
        my $hseq2 = substr($current_seq, $cropstop);
        if ( !($pseq) || !($hseq1) || !($hseq2)){
            err("ERROR: problem with mid cropping $current_contig\n");
        }
        print_seq($OA, "$current_contig\_HAPLOTIG_C1", $hseq1);
        print_seq($OA, "$current_contig\_HAPLOTIG_C2", $hseq2);
        print_seq($OG, $current_contig, $pseq);
        
        # left-crop
    } elsif ($cropstart == 0){
        my $pseq = substr($current_seq, $cropstart, $cropstop);
        my $hseq = substr($current_seq, $cropstop);
        if (!($pseq) || !($hseq)){
            err("ERROR: problem left cropping $current_contig\n");
        }
        print_seq($OA, "$current_contig\_HAPLOTIG_C1", $hseq);
        print_seq($OG, $current_contig, $pseq);
        
        # right-crop
    } else {
        my $pseq = substr($current_seq, $cropstart);
        my $hseq = substr($current_seq, 0, $cropstart);
        if (!($pseq) || !($hseq)){
            err("ERROR: problem right cropping $current_contig\n");
        }
        print_seq($OA, "$current_contig\_HAPLOTIG_C1", $hseq);
        print_seq($OG, $current_contig, $pseq);
    }
}

sub print_seq{
    my $FH = $_[0];
    my $id = $_[1];
    my $seq = $_[2];
    
    print $FH ">$id";
    
    for (my $i=0; $i < (length($seq)-1); $i += 100){
        print $FH "\n";
        print $FH substr($seq, $i, 100);
    }
    
    print $FH "\n";
}


sub err {
    print STDERR @_;
    exit(1);
}



