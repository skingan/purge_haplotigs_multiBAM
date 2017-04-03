#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);


#---PARAMS---

my $threads = 4;
my $maxmatch_cutoff = 175;
my $bestmatch_cutoff = 85;

my $stats_csv;
my $genome_fasta;
my $unknown_only;

my $script_dir = "$Bin/../scripts";
my $temp_dir = "tmp_reassign_contigs";
my $minced_dir = "tmp_reassign_contigs/minced_genome";


#---HELP MSG---

my $usage = "
Usage:
zz2_assign_contigs.pl  -s  stats.csv  -g  genome.fasta  [ -p -m -a ]

REQUIRED:
-s      The stats .csv file output from zz1_analyze_genecov.pl
-g      Genome .fasta (also need the samtools index genome.fasta.fai file)

OPTIONAL:
-p      Threads to use for the blastn search. DEFAULT = $threads
-m      Maxmatch cutoff percentage. Used to determine if a contig is a 
        repetitive sequence. DEFAULT = $maxmatch_cutoff
-a      Bestmatch cutoff percentage. Used to determine if a contig is a
        haplotig. DEFAULT = $bestmatch_cutoff
-u      Produce dotplots for unknown contigs only (dont produce them
        for the auto-assigned contigs)
";

#---PARSE ARGS---

GetOptions(
    "stats=s" => \$stats_csv,
    "genome=s" => \$genome_fasta,
    "proc=i" => \$threads,
    "maxmatch=i" => \$maxmatch_cutoff,
    "allmatch=i" => \$bestmatch_cutoff,
    "unknown" => \$unknown_only
) or die $usage;

if ( !($stats_csv) || !($genome_fasta) ){
    err($usage);
}


#---CHECK DEPENDENCIES---

my $check_die = 0;

my $blastn = `blastn -version`;
if ($blastn eq ""){
    print STDERR "ERROR: blastn doesn't appear to be installed\n";
    $check_die = 1;
}

# TODO check scripts and mummer

check_files($stats_csv, $genome_fasta);

if ($check_die){
    err("ERROR: One or more errors encountered, exiting\n");
}


#---PIPELINE---

# BLASTN DB
if (!(-d "$temp_dir/blstdb")){
    print STDERR "INFO: blastn db not found, creating\n";
    mkdir "$temp_dir/blstdb";
    runcmd("makeblastdb -in $genome_fasta -dbtype nucl -out $temp_dir/blstdb/$genome_fasta\n");
} else {
    print STDERR "INFO: blastn db found, skipping\n";
}

# list of suspects and junk
if (!(-s "$temp_dir/suspects.list")){
    print STDERR "INFO: suspects.list not found, creating\n";
    my @suspects;
    my @junk;
    open(my $CSV, $stats_csv) or err("ERROR: Failed to open $stats_csv for reading\n");
    while(<$CSV>){
        next if ($_ =~ /^#/);
        if ($_ =~ /([a-zA-Z0-9_-]+),[sS],/){
            push @suspects, $1;
        } elsif ($_ =~ /([a-zA-Z0-9_-]+),[jJ],/){
            push @junk, $1;
        } elsif ($_ !~ /([a-zA-Z0-9_-]+),/){
            err("ERROR: Failed to get seq id from entry $_, illegal characters in seq name? Seq id should contain only a-z, A-Z, 0-9, _, -\n");
        }
    }
    open my $SL, ">", "$temp_dir/suspects.list" or err("ERROR: Failed to open $temp_dir/suspects.list for writing\n");
    foreach(@suspects){
        print $SL "$_\n";
    }
    close($SL);
    open my $JL, ">", "$temp_dir/junk.list" or err("ERROR: Failed to open $temp_dir/junk.list for writing\n");
    foreach(@junk){
        print $JL "$_\n";
    }
    close($JL);
} else {
    print STDERR "INFO: suspects.list found, skipping\n";
}


# do the suspects v all blastn
if (!(-s "$temp_dir/suspects.blastn.gz")){
    print STDERR "INFO: suspects.blastn.gz not found, creating\n";
    # get the suspect contig sequences
    if (!(-s "$temp_dir/suspects.fasta")){
        print STDERR "INFO: suspects.fasta not found, needed for blastn, creating\n";
        runcmd("$script_dir/returnseq.pl -f $genome_fasta -l $temp_dir/suspects.list > $temp_dir/suspects.fasta\n");
    } else {
        print STDERR "INFO: suspects.fasta found, skipping\n";
    }
    runcmd("blastn -query $temp_dir/suspects.fasta -db $temp_dir/blstdb/$genome_fasta -outfmt 6 -num_alignments 3 -evalue 0.0000001 -num_threads $threads |  awk ' \$1 != \$2 && \$4 > 500 { print } ' | gzip - > $temp_dir/suspects.blastn.gz\n");
    unlink "$temp_dir/suspects.fasta";
} else {
    print STDERR "INFO: suspects.blastn.gz found, skipping\n";
}

# get the chopped up suspects for speed
if (!(-d $minced_dir)){
    print STDERR "INFO: $minced_dir not found, creating and mincing genome\n";
    mkdir $minced_dir;
    runcmd("$script_dir/mince_genome.pl -f $genome_fasta -o $minced_dir\n");
}


# analyse the blastn output, gen dotplots for reciprocal best hits
print STDERR "INFO: Starting analysis and assigning contigs, this may take a while, check $temp_dir/analyse_blastn.log and suspect_contig_reassign.tsv for progress\n";
if ($unknown_only){
    runcmd("$script_dir/analyze_blastn_output.pl -f $genome_fasta.fai -d $minced_dir -b $temp_dir/suspects.blastn.gz -t suspect_contig_reassign.tsv -u 1> $temp_dir/analyse_blastn.log 2>&1 \n");
} else {
    runcmd("$script_dir/analyze_blastn_output.pl -f $genome_fasta.fai -d $minced_dir -b $temp_dir/suspects.blastn.gz -t suspect_contig_reassign.tsv 1> $temp_dir/analyse_blastn.log 2>&1 \n");
}




#---SUBROUTINES---

sub err {
    print STDERR @_;
    exit(1);
}

sub runcmd {
  print STDERR "RUNNING: @_";
  system(@_) == 0 or err("ERROR: Failed to run command: ", @_);
}

sub check_files {
    foreach(@_){
        if (!(-e $_)){
            print STDERR "ERROR: No such file, or empty file for $_\n";
            $check_die = 1;
        }
    }
}

