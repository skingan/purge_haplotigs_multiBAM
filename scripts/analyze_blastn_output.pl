#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use threads;
use Thread::Semaphore;


# params
my $fasta;
my $blast;
my $tblout;
my $seq_dir;
my $unknown_only;
my $max_match_cutoff = 250;
my $align_match_cutoff = 75;
my $threads = 4;
my $no_call_cutoff = 20;

my $usage = "
Usage:
analyze_blastn_output.pl  -f genome.fasta  -b blastn_out.fmt6.gz  -t CONTIG_REASSIGN.tsv  [ -p $threads -m $max_match_cutoff  -a $align_match_cutoff ]

-f      genome.fasta file, needs to be indexed with samtools faidx
-d      suspect fasta seq directory, output from zz2_assign_contigs.pl
-b      the output from zz1_ pipeline step (blastn, -outfmt 6, optionally gzipped)
-t      output .tsv file which then needs to be edited and fed into next step

OPTIONAL:
-p      Threads, default = $threads.
-m      maxmatch cutoff - for finding collapsed repeat contigs. -m $max_match_cutoff means that if a contig 
            matches $max_match_cutoff % or more to another contig then it's a collased repeat/assembly junk.
-a      best alignment cutoff - for identifying haplotigs. -a $align_match_cutoff means that if a contig is not
            a collapsed repeat, and maps $align_match_cutoff % or more to another contig, then it's a haplotig.
-u      Produce dotplots for unknown suspects only (don't make dotplots for auto-assigned contigs);
";


GetOptions(
    "fastafai=s" => \$fasta,
    "dir=s" => \$seq_dir,
    "blast=s" => \$blast,
    "table=s" => \$tblout, 
    "unknown" => \$unknown_only,
    "procs=i" => \$threads,
    "maxmatch=s" => \$max_match_cutoff,
    "alignmatch=s" => \$align_match_cutoff
) or err($usage);

if (!($fasta) || !($blast) || !($tblout) || !($seq_dir)){
    err($usage);
}

# genome index
if (!(-s $fasta)){
    err("$fasta needs to be the samtools faidx index file for the genome, i.e. genome.fasta.fai\nERROR: cant find $fasta\n");
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
if (!(-d $dotcall) && !($unknown_only)){
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


# set up for multithreading
my $available_threads = Thread::Semaphore->new($threads);
my $writing_to_out = Thread::Semaphore->new(1);
my $jobnumber = 0;

# open files
open($FAI, $fasta) or err("failed to open $fasta for reading\n");
if ($blast =~ /\.gz$/){
    open($BLST, "gunzip -c $blast |") or err("failed to open gunzip pipe from $blast\n");
} else {
    open($BLST, $blast) or err("failed to open file $blast for reading\n");
}
open($OUT, ">$tblout") or err("failed to open $tblout for writing\n");

# read in the contig lengths
while (<$FAI>){
    my @line = split(/\s+/, $_);
    $line[0] =~ s/\|.+//;
    $length{$line[0]} = $line[1];
}
close($FAI);

# read in blast output file
while (<$BLST>){
    my @line = split(/\s+/, $_);
    $line[0] =~ s/\|.+//;
    $line[1] =~ s/\|.+//;
    if (!($hits{$line[0]}{1})){
        $hits{$line[0]}{1} = $line[1] if ($line[0] ne $line[1]); # shouldn't be self hits but just in case
    } elsif (!($hits{$line[0]}{2}) && ($hits{$line[0]}{1} ne $line[1])){
        $hits{$line[0]}{2} = $line[1] if ($line[0] ne $line[1]);
    }
}
close($BLST);

# characterise hits
foreach my $contig (keys(%hits)){
    # check for length
    if (!($length{$contig})){
        err("ERROR: missing length in genome .fasta.fai index file for $contig\n");
    }
    
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
    
    # ignore if contig is larger that both its hit seqs (or one best hit if it only has one)
    elsif ($hits{$contig}{2}) {
        if ( $length{$contig} > ($length{$hits{$contig}{1}} + $length{$hits{$contig}{2}}) ){
            $purged{$contig} = 1;
        }
    }
    elsif ($length{$contig} > $length{$hits{$contig}{1}}) {
        $purged{$contig} = 1;
    }
}

if ($unknown_only){
    print STDERR "INFO: skipping generation of dotplots for auto-assigned contigs\n";
}

print $OUT "#for cropping, specify the section to KEEP\n";
print $OUT "#example_crop1\tmatch1\tmatch2\t50.00\t50.00\tc\tstart\t500000\n";
print $OUT "#example_crop2\tmatch1\tmatch2\t50.00\t50.00\tc\t500000\tend\n";
print $OUT "#example_crop3\tmatch1\tmatch2\t50.00\t50.00\tc\t250000\t750000\n";
print $OUT "#note for reassign key: everything other than 'c', 'r', and 'h' will be ignored\n";
print $OUT "#\n";
print $OUT "#SUSPECT_CONTIG\tTOP_MATCH\tSECOND_MATCH\tMaxMatchCov\tBestMatchCov\tREASSGIN_KEY\tCROP_START\tCROP_END\n";
print $OUT "#Reciprocal_best_hits\n";

foreach my $contig (keys(%RBHs)){
    $available_threads->down(1);
    threads->create(\&RBH, $jobnumber, $contig);
    $jobnumber++;
}

print $OUT "#Everything_else\n";
foreach my $contig (sort(keys(%hits))){
    $available_threads->down(1);
    threads->create(\&EE, $jobnumber, $contig);
    $jobnumber++;
}

# wait on remaining jobs
while(){
    if (threads->list()){
        sleep 0.5;
    } else {
        last;
    }
}

close($OUT);

exit(0);


#---SUBROUTINES---

sub RBH {
    my $job = $_[0];
    my $contig = $_[1];
    
    if (!($purged{$RBHs{$contig}})){
        # we assume the smaller contig is the haplotig
        runcmd("cat $seq_dir/$hits{$RBHs{$contig}}{1}.fasta > $temp/$job.ref.fasta\n");
        if ($hits{$RBHs{$contig}}{2}){
            runcmd("cat $seq_dir/$hits{$RBHs{$contig}}{2}.fasta >> $temp/$job.ref.fasta\n");
        }
        my $assign = guess_assignment($RBHs{$contig}, $job);
        $writing_to_out->down(1);
        print $OUT "$RBHs{$contig}\t$contig\t$hits{$RBHs{$contig}}{2}\t$assign\n";
        $writing_to_out->up(1);
    }
    
    # exit
    $available_threads->up(1);
    threads->detach();
}

sub EE {
    my $job = $_[0];
    my $contig = $_[1];
    
    if (!($purged{$contig}) && ($hits{$contig}{1})){
        runcmd("cat $seq_dir/$hits{$contig}{1}.fasta > $temp/$job.ref.fasta\n");
        if ($hits{$contig}{2}){
            runcmd("cat $seq_dir/$hits{$contig}{2}.fasta >> $temp/$job.ref.fasta\n");
            my $assign = guess_assignment($contig, $job);
            $writing_to_out->down(1);
            print $OUT "$contig\t$hits{$contig}{1}\t$hits{$contig}{2}\t$assign\n";
            $writing_to_out->up(1);
        }
    }
    # exit
    $available_threads->up(1);
    threads->detach();
}

sub guess_assignment{
    my $query = $_[0];
    my $job = $_[1];
    
    my @LOG = ();
    
    my $ref = "$temp/$job.ref.fasta";
    # get seq
    runcmd("cat $seq_dir/$query.fasta > $temp/$job.tmp_query.fasta\n");
    
    # nucmer
    runcmd("$MDPbin/MDP_nucmer --maxmatch -p $temp/$job.tmp $temp/$job.tmp_query.fasta $ref\n");
    
    # deltas
    runcmd("$MDPbin/MDP_delta-filter -m $temp/$job.tmp.delta > $temp/$job.tmp.m.delta\n");
    runcmd("$MDPbin/MDP_delta-filter -r $temp/$job.tmp.delta > $temp/$job.tmp.r.delta\n");
    
    # all-match coverage
    push @LOG, "$MDPbin/MDP_show-coords -b -c $temp/$job.tmp.m.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'\n";
    my $maxmatch = `$MDPbin/MDP_show-coords -b -c $temp/$job.tmp.m.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'`;
    $maxmatch =~ s/\s//g;
    push @LOG, "MAXMATCH coverage = $maxmatch\n";
    
    # best-align coverage
    push @LOG, "$MDPbin/MDP_show-coords -b -c $temp/$job.tmp.r.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'\n";
    my $alignmatch = `$MDPbin/MDP_show-coords -b -c $temp/$job.tmp.r.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'`;
    $alignmatch =~ s/\s//g;
    push @LOG, "BESTMATCH coverage = $alignmatch\n";
    
    my $a = "$maxmatch\t$alignmatch\t";
    
    # guess the assignment and make dotplots
    if ( ($maxmatch >= $max_match_cutoff) && ($alignmatch >= $align_match_cutoff) ) {
        $purged{$query}=1;
        $a .= "r";
        push @LOG, "\n### $query = repeat/assembly junk\n\n";
        runcmd("$MDPbin/MDP_mummerplot --fat -p $dotcall/$query $temp/$job.tmp.m.delta") if (!($unknown_only));
        
    } elsif ($alignmatch >= $align_match_cutoff){
        $purged{$query}=1;
        $a .= "h";
        push @LOG, "\n### $query = haplotig\n\n";
        runcmd("$MDPbin/MDP_mummerplot --fat -p $dotcall/$query $temp/$job.tmp.m.delta") if (!($unknown_only));
        
    } elsif ($maxmatch < $no_call_cutoff){
        $a .= "?";
        push @LOG, "\n### $query = no call, insufficient seq homology to reference hits\n\n";
    } else {
        $purged{$query}=1;
        $a .= "?";
        push @LOG, "\n### $query = unsure, use your mk-I eyeballs\n\n";
        runcmd("$MDPbin/MDP_mummerplot --fat -p $dotunk/$query $temp/$job.tmp.m.delta\n");
    }
    
    # clean-up
    my @files_to_clean_up = (
        "$dotcall/$query.filter", "$dotcall/$query.fplot", "$dotcall/$query.rplot", "$dotcall/$query.gp",
        "$dotunk/$query.filter", "$dotunk/$query.fplot", "$dotunk/$query.rplot", "$dotunk/$query.gp",
        "$temp/$job.tmp_query.fasta", "$temp/$job.tmp.delta", "$temp/$job.tmp.m.delta", "$temp/$job.tmp.r.delta", "$temp/$job.ref.fasta"
    );
    foreach my $file (@files_to_clean_up){
        if (-s $file){
            push @LOG, "INFO: cleaning up $file\n";
            unlink $file;
        }
    }
    $writing_to_out->down(1);
    print STDERR @LOG;
    $writing_to_out->up(1);
    
    return($a);
}

sub err {
    print STDERR @_;
    exit(1);
}


sub runcmd {
  print STDERR "RUNNING: @_";
  system(@_) == 0 or err("ERROR: Failed to run command: ", @_);
}
