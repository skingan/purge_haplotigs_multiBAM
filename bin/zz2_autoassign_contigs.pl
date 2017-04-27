#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use Time::Piece;


#---PARAMS---

my $threads = 4;
my $maxmatch_cutoff = 250;
my $bestmatch_cutoff = 75;
my $low_cutoff = 40;
my $passes = 3;

my $stats_csv;
my $genome_fasta;
my $unknown_only;
my $current_pass = 1;

my $script_dir = "$Bin/../scripts";
my $temp_dir = "tmp_reassign_contigs";
my $minced_dir = "tmp_reassign_contigs/minced_genome";
my $outfile = "suspect_contig_reassign.tsv";
my $curated = "curated";
my @haplotig_files;
my @artefact_files;
my @reassignment_logs;
my %ignore;

#---HELP MSG---

my $usage = "
Usage:
zz2_assign_contigs.pl  -s  stats.csv  -g  genome.fasta  [ -o -t -p -c -m -a ]

REQUIRED:
-s      The stats .csv file output from zz1_analyze_genecov.pl
-g      Genome .fasta (also need the samtools index genome.fasta.fai file)

OPTIONAL:
-o      Output file name for the reassignment .tsv file. 
        DEFAULT = $outfile
-t      Threads to use for the blastn search and mummer alignments.
        DEFAULT = $threads
-p      Max number of passes to perform, DEFAULT = $passes. More than one purging
        pass is usually needed due to overlapping haplotigs and repeat contigs.
-c      Prefix for the curated assembly, DEFAULT = $curated.

-m      Maxmatch cutoff percentage. Used to determine if a contig is a 
        repetitive sequence. DEFAULT = $maxmatch_cutoff
-a      Bestmatch cutoff percentage. Used to determine if a contig is a
        haplotig. DEFAULT = $bestmatch_cutoff

-u      Produce dotplots for unknown contigs only. DEFAULT = produce
        dotplots for both assigned and unassigned.
";


msg("Starting pipeline using command:\n$0 @ARGV\n");


#---PARSE ARGS---

GetOptions(
    "stats=s" => \$stats_csv,
    "genome=s" => \$genome_fasta,
    "threads=i" => \$threads,
    "passes=i" => \$passes,
    "curated=s" => \$curated,
    "outfile=s" => \$outfile,
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
    print STDERR "ERROR: blastn doesn't appear to be installed";
    $check_die = 1;
}

my $gnuparallel = `parallel --help`;
if ($gnuparallel eq ""){
    $gnuparallel = 0;
} else {
    $gnuparallel = 1;
}

# TODO check scripts and mummer

check_files($stats_csv, $genome_fasta);

if ($check_die){
    err("ERROR: One or more errors encountered, exiting\n");
}

if (!(-d $temp_dir)){
    mkdir $temp_dir;
}

suspects_and_junk();

mince_genome();


#---PIPELINE---

PURGE: while(1){
    
    msg("\n\n#####\n\nPerforming purging pass $current_pass\n\n#####\n\n");
    
    refine_suspects();
    
    blastn_db();
    
    suspects_blastn();
    
    analyse_blastn();
    
    reassign_contigs();
    
    my $prev = $genome_fasta;
    $genome_fasta = "pass_$current_pass.fasta";
    
    my @sp = stat $prev;
    my @sc = stat "pass_$current_pass.fasta";
    
    if ( ($current_pass == $passes) || ($sp[7] == $sc[7]) ){
        last PURGE;
    } else {
        reset_purging();
    }
    
    $current_pass++;
}

rename_output();


# done
msg("Contig assignment completed successfully!\n");


#---SUBROUTINES---

sub suspects_and_junk {
    # list of suspects and junk
    if (!(-s "$temp_dir/suspects.list")){
        msg("suspects.list not found, creating");
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
        msg("suspects.list found, skipping");
    }
}

sub mince_genome {
    # get the chopped up suspects for speed
    if (!(-d $minced_dir)){
        msg("$minced_dir not found, creating and mincing genome");
        mkdir $minced_dir;
        runcmd("$script_dir/mince_genome.pl -f $genome_fasta -o $minced_dir\n");
    }
}

#---MULTI PASS SUBS---

sub blastn_db {
    # BLASTN DB
    if (!(-d "$temp_dir/blstdb")){
        msg("blastn db not found, creating");
        mkdir "$temp_dir/blstdb";
        runcmd("makeblastdb -in $genome_fasta -dbtype nucl -out $temp_dir/blstdb/$genome_fasta\n");
        print STDERR "\n\n";
    } else {
        msg("blastn db found, skipping");
    }
}

sub refine_suspects {
    if ($current_pass > 1){
        # get the contigs flagged for 'no_reassignment' due to poor seq homology
        open(my $OF, "$outfile") or err("Failed to open $outfile for reading");
        while(<$OF>){
            if ($_ =~ /NO_REASSIGNMENT/){
                my @line = split(/\s+/, $_);
                $ignore{$line[0]}=1;
            }
        }
        close($OF);
        # remake the suspects.list to remove 'no_reassignment' contigs to improve speed with subsequent passes
        my @clean_suspects;
        open(my $SL, "$temp_dir/suspects.list") or err("Failed to open $temp_dir/suspects.list for reading");
        while(<$SL>){
            $_ =~ s/\s//g;
            if (!($ignore{$_})){
                push @clean_suspects, $_;
            }
        }
        close($SL);
        open($SL, ">", "$temp_dir/suspects.list") or err("Failed to open $temp_dir/suspects.list for writing");
        (print $SL "$_\n") foreach(@clean_suspects);
        close($SL);
        undef @clean_suspects;
        undef %ignore;
    }
}

sub suspects_blastn {
    # do the suspects v all blastn
    if (!(-s "$temp_dir/suspects.blastn.gz")){
        msg("suspects.blastn.gz not found, creating");
        # get the suspect contig sequences
        if (!(-s "$temp_dir/suspects.fasta")){
            msg("suspects.fasta not found, needed for blastn, creating");
            runcmd("$script_dir/returnseq.pl -f $genome_fasta -l $temp_dir/suspects.list > $temp_dir/suspects.fasta");
        } else {
            msg("suspects.fasta found, skipping");
        }
        if ($gnuparallel){
            runcmd("cat $temp_dir/suspects.fasta | parallel --no-notice -j $threads --block 100k --recstart '>' --pipe blastn -db $temp_dir/blstdb/$genome_fasta -outfmt 6 -evalue 0.000000000001 -max_target_seqs 3 -max_hsps 1000 -word_size 28 -culling_limit 10 -query - | awk ' \$1 != \$2 && \$4 > 500 { print } ' | gzip - > $temp_dir/suspects.blastn.gz");
        } else {
            runcmd("blastn -query $temp_dir/suspects.fasta -db $temp_dir/blstdb/$genome_fasta -outfmt 6 -evalue 0.000000000001 -num_threads $threads -max_target_seqs 3 -max_hsps 1000 -word_size 28 -culling_limit 10 |  awk ' \$1 != \$2 && \$4 > 500 { print } ' | gzip - > $temp_dir/suspects.blastn.gz");
        }
        unlink "$temp_dir/suspects.fasta";
    } else {
        msg("suspects.blastn.gz found, skipping");
    }
}

sub analyse_blastn {
    # analyse the blastn output, gen dotplots for reciprocal best hits
    msg("Starting analysis and assigning contigs, this may take a while, check $temp_dir/analyse_blastn.log and $outfile for progress");
    if ($unknown_only){
        runcmd("$script_dir/analyze_blastn_output.pl -f $genome_fasta.fai -d $minced_dir -b $temp_dir/suspects.blastn.gz -t $outfile -m $maxmatch_cutoff -a $bestmatch_cutoff -u -p $threads 1> $temp_dir/analyse_blastn.log 2>&1");
    } else {
        runcmd("$script_dir/analyze_blastn_output.pl -f $genome_fasta.fai -d $minced_dir -b $temp_dir/suspects.blastn.gz -t $outfile -m $maxmatch_cutoff -a $bestmatch_cutoff -p $threads 1> $temp_dir/analyse_blastn.log 2>&1");
    }
}

sub reassign_contigs {
    runcmd("$Bin/zz3_reassign_contigs.pl -t $outfile -g $genome_fasta -o pass_$current_pass");
    runcmd("samtools faidx pass_$current_pass.fasta");
    push @haplotig_files, "pass_$current_pass.haplotigs.fasta";
    push @artefact_files, "pass_$current_pass.artefacts.fasta";
    push @reassignment_logs, "pass_$current_pass.reassignments.tsv";
}

sub reset_purging {
    runcmd("cp $outfile pass_$current_pass.$outfile");
    foreach my $file ("$temp_dir/suspects.blastn.gz", "$temp_dir/analyse_blastn.log"){
        msg("Cleaning up temp file $file");
        unlink $file or err("Failed to clean up temp file");
    }
    foreach my $dir ("$temp_dir/blstdb", "dotplots_unknowns"){
        if (-d $dir){
            msg("Cleaning up directory $dir");
            runcmd("rm -rf $dir");
        }
    }
}

sub rename_output {
    runcmd("mv pass_$current_pass.fasta $curated.fasta");
    runcmd("mv pass_$current_pass.fasta.fai $curated.fasta.fai");
    my $hf = "$curated.haplotigs.fasta";
    my $af = "$curated.artefacts.fasta";
    my $al = "$curated.reassignments.tsv";
    
    if (-s $hf){
        unlink $hf;
    }
    foreach my $file (@haplotig_files){
        runcmd("cat $file >> $hf");
    }
    if (-s $af){
        unlink $af;
    }
    foreach my $file (@artefact_files){
        runcmd("cat $file >> $af");
    }
    if (-s $al){
        unlink $al;
    }
    foreach my $file (@reassignment_logs){
        runcmd("cat $file >> $al");
    }
    
    # get the contig association paths for the reassigned contigs
    runcmd("$script_dir/get_reassignment_paths.pl -i $al -o $curated.reassignment_paths.log");
}

#---UTILITY SUBROUTINES---

sub msg {
    my $t = localtime;
    my $line = $t->dmy . " " . $t->hms . " INFO: @_\n";
    print STDERR $line;
}

sub err {
    msg(@_);
    exit(1);
}

sub runcmd {
  msg("RUNNING: @_");
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

