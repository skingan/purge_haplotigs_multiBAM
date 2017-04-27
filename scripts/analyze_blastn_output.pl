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
my $no_call_cutoff = 40;

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
    
    # ignore if contig is larger than both its hit seqs (or one best hit if it only has one)
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

foreach my $contig (sort(keys(%hits))){
    $available_threads->down(1);
    threads->create(\&ASSIGN, $jobnumber, $contig);
    $jobnumber++;
}


# wait on remaining jobs
sleep 1 while threads->list();


close($OUT);

exit(0);


#---SUBROUTINES---

sub ASSIGN {
    my $job = $_[0];
    my $contig = $_[1];
    my $cmd;
    
    my $LOG;
    
    if (!($purged{$contig}) && ($hits{$contig}{1})){
        
        $cmd = "cat $seq_dir/$hits{$contig}{1}.fasta > $temp/$job.ref.fasta\n";
        $LOG .= $cmd;
        runcmd($cmd);
        
        if ($hits{$contig}{2}){
            
            $cmd = "cat $seq_dir/$hits{$contig}{2}.fasta >> $temp/$job.ref.fasta\n";
            $LOG .= $cmd;
            runcmd($cmd);
        }
        
        my ($assign, $S_LOG) = guess_assignment($contig, $job);
        
        $LOG .= $S_LOG;
        
        # print output
        $writing_to_out->down(1);
        print $OUT "$contig\t$hits{$contig}{1}\t$hits{$contig}{2}\t$assign\n";
        # print log
        print STDERR "$LOG\n###\n\n";
        $writing_to_out->up(1);
    }
    # exit
    $available_threads->up(1);
    threads->detach();
}

sub guess_assignment{
    my $query = $_[0];
    my $job = $_[1];
    my $LOG;
    my $cmd;
    
    my $ref = "$temp/$job.ref.fasta";
    # get seq
    $cmd = "cat $seq_dir/$query.fasta > $temp/$job.tmp_query.fasta\n";
    $LOG .= $cmd;
    runcmd($cmd);
    
    # nucmer
    $cmd = "$MDPbin/MDP_nucmer --maxmatch -p $temp/$job.tmp $temp/$job.tmp_query.fasta $ref 2>&1\n";
    $LOG .= $cmd;
    $LOG .= `$cmd`;
    
    # deltas
    $cmd = "$MDPbin/MDP_delta-filter -m $temp/$job.tmp.delta > $temp/$job.tmp.m.delta\n";
    $LOG .= $cmd;
    runcmd($cmd);
    
    $cmd = "$MDPbin/MDP_delta-filter -r $temp/$job.tmp.delta > $temp/$job.tmp.r.delta\n";
    $LOG .= $cmd;
    runcmd($cmd);
    
    # all-match coverage
    $cmd = "$MDPbin/MDP_show-coords -b -c $temp/$job.tmp.m.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'\n";
    $LOG .= $cmd;
    my $maxmatch = `$cmd`;
    $maxmatch =~ s/\s//g;
    $LOG .= "MAXMATCH coverage = $maxmatch\n";
    
    # best-align coverage
    $cmd = "$MDPbin/MDP_show-coords -b -c $temp/$job.tmp.r.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'\n";
    $LOG .= $cmd;
    my $alignmatch = `$cmd`;
    $alignmatch =~ s/\s//g;
    $LOG .= "BESTMATCH coverage = $alignmatch\n";
    
    my $a = "$maxmatch\t$alignmatch\t";
    
    # guess the assignment and make dotplots
    if ( ($maxmatch >= $max_match_cutoff) && ($alignmatch >= $align_match_cutoff) ) {
        $purged{$query}=1;
        $a .= "r";
        $LOG .= "\n### $query = repeat/assembly junk\n\n";
        
        if (!($unknown_only)){
            $cmd = "$MDPbin/MDP_mummerplot --fat -p $dotcall/$query $temp/$job.tmp.m.delta 2>&1 > /dev/null\n";
            $LOG .= $cmd;
            runcmd($cmd);
        }
        
    } elsif ($alignmatch >= $align_match_cutoff){
        $purged{$query}=1;
        $a .= "h";
        $LOG .= "\n### $query = haplotig\n\n";
        
        if (!($unknown_only)){
            $cmd = "$MDPbin/MDP_mummerplot --fat -p $dotcall/$query $temp/$job.tmp.m.delta 2>&1 > /dev/null\n";
            $LOG .= $cmd;
            runcmd($cmd);
        }
        
    } elsif ($alignmatch < $no_call_cutoff){
        $a .= "NO_REASSIGNMENT";
        $LOG .= "\n### $query = no call, insufficient seq homology to reference hits\n\n";
    } else {
        $purged{$query}=1;
        $a .= "?";
        $LOG .= "\n### $query = unsure, use your mk-I eyeballs\n\n";
        
        $cmd = "$MDPbin/MDP_mummerplot --fat -p $dotunk/$query $temp/$job.tmp.m.delta 2>&1 > /dev/null\n";
        $LOG .= $cmd;
        runcmd($cmd);
    }
    
    # clean-up
    my @files_to_clean_up = (
        "$dotcall/$query.filter", "$dotcall/$query.fplot", "$dotcall/$query.rplot", "$dotcall/$query.gp",
        "$dotunk/$query.filter", "$dotunk/$query.fplot", "$dotunk/$query.rplot", "$dotunk/$query.gp",
        "$temp/$job.tmp_query.fasta", "$temp/$job.tmp.delta", "$temp/$job.tmp.m.delta", "$temp/$job.tmp.r.delta", $ref
    );
    
    foreach my $file (@files_to_clean_up){
        if (-s $file){
            $LOG .= "INFO: cleaning up $file\n";
            unlink $file;
        }
    }
    
    return($a, $LOG);
}

sub err {
    print STDERR @_;
    exit(1);
}


sub runcmd {
  system(@_) == 0 or err("ERROR: Failed to run command: ", @_);
}
