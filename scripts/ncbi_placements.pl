#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use threads;
use threads::shared;
use Thread::Semaphore;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use PipeUtils;
use List::Util qw(min max);


# Adaptation of nucmer2ncbiPlacement.py from https://github.com/skingan/NCBI_DiploidAssembly/blob/master/nucmer2ncbiPlacement.py
# The adaptation is necessary if the assembly has been curated with purge_haplotigs as the FALCON IDs won't necessarily match up.


# PIPELINE:
    # blastn haplotigs against primary contigs and get best hit for each
    # nucmer, delta-filter, show-coords each haplotig against it's blastn hit
    # get the alignment information from the coords file and write to the output placements file


#---SET UP LOGGING---
our $LOG;
my $TMP = "temp_ncbi_placements";
(-d $TMP) or mkdir $TMP;
open $LOG, ">", "$TMP/purge_haplotigs_ncbiplace.log" or err("failed to open log file for writing");



#---GLOBAL VARIABLES---
my $primary;
my $haplotigs;
my $out = "ncbi_placements.tsv";
my $threads = 4;
my $coverage = 50;

my %output :shared;     # $output{$htig}{len} = h-tig length
                        #               {top} = top hit
                        #               {rlen} = ref alignments running total (negative = rev alignments)
                        #               {hlen} = htig alignments running total
                        #               {hs} = htig start
                        #               {he} = htig end
                        #               {rs} = ref start
                        #               {re} = ref end


my $usage = "
USAGE:
purge_haplotigs  ncbiplace  -p primary_contigs.fasta  -h haplotigs.fasta  [ -o outfile  -t threads  -c coverage]

REQUIRED:
-p / -primary       Primary contigs fasta file
-h / -haplotigs     Haplotigs fasta file

OPTIONAL:
-o / -out           Out file name (DEFAULT = $out)
-t / -threads       Threads for blastn and nucmer steps (DEFAULT = $threads)
-c / -coverage      Coverage cutoff percentage for pairing contigs (DEFAULT = $coverage %)

IMPORTANT NOTES:  
    Will only work with nucmer v3.9 or higher

    The placement file generated is currently untested and this script is still under development.
    This script is intended to be an adaptation of nucmer2ncbiPlacement.py from:
    https://github.com/skingan/NCBI_DiploidAssembly/blob/master/nucmer2ncbiPlacement.py
    The adaptation is necessary if the assembly has been curated with Purge Haplotigs as the 
    curated primary contigs and haplotigs wont necessarily be paired by the FALCON Unzip IDs.

";


#---PREFLIGHT CHECKS---
check_programs("blastn", "samtools", "parallel", "nucmer", "delta-filter", "show-coords");


#---PARSE ARGS---
my $args = "@ARGV";
GetOptions (
    "primary=s" => \$primary,
    "haplotigs=s" => \$haplotigs,
    "out=s" => \$out,
    "threads=i" => \$threads,
    "coverage=i" => \$coverage
) or die $usage;


# check args and files
($primary) && ($haplotigs) || die $usage;

check_files($primary, $haplotigs);



#---THREADS---

my $available_threads = Thread::Semaphore->new($threads);
my $writing_to_out = Thread::Semaphore->new(1);



# print settings
msg("

PARAMETERS:
Primary contigs FASTA       $primary
Haplotigs FASTA             $haplotigs
Out FASTA                   $out
Threads                     $threads
Coverage align cutoff       $coverage %

RUNNING USING COMMAND:
purge_haplotigs ncbiplace $args
");


# make the temp dir
mkdir $TMP if (!(-d $TMP));

# get fasta indexes ready and get the haplotig IDs and lengths
((-s "$_.fai") || runcmd({ command => "samtools faidx $_ 2> $TMP/samtools_faidx.stderr", logfile => "$TMP/samtools_faidx.stderr" })) for ($haplotigs, $primary);

# read in the contig lengths
open my $HFI, "<", "$haplotigs.fai" or err("failed to open file $haplotigs.fai");
while(<$HFI>){
    my @l = split(/\s+/, $_);
    my %l : shared;
    $output{$l[0]} = \%l;
    $output{$l[0]}{len} = $l[1];
}
close $HFI;



#---BLASTN TO GET HITS FOR ALL HAPLOTIGS---
if (!(-s "$TMP/$primary.outfmt6")){
    runcmd({ command => "cat $haplotigs | parallel --no-notice -j $threads --block 100k --recstart '>' --pipe blastn -subject $primary -outfmt 6 -evalue 0.00000001 -max_target_seqs 1 -max_hsps 1000 -word_size 28 -culling_limit 10 -query - > $TMP/$primary.outfmt6.tmp 2> $TMP/blast.stderr  &&  mv  $TMP/$primary.outfmt6.tmp $TMP/$primary.outfmt6", logfile => "$TMP/blast.stderr" });
} else {
    msg("found $TMP/$primary.outfmt6, SKIPPING blastn step");
}

# get an index of the hits
open my $BLST, "<", "$TMP/$primary.outfmt6" or err("failed to open $TMP/$primary.outfmt6 for reading");
while(<$BLST>){
    my @l = split(/\s+/, $_);
    ($output{$l[0]}{top}) || ($output{$l[0]}{top} = $l[1]);
}
close $BLST;

# open outfile for writing
open my $OUT, ">", $out or err("failed to open $out for writing");
# print the output header
print $OUT "#alt_asm_name\tprim_asm_name\talt_scaf_name\tparent_type\tparent_name\tori\talt_scaf_start\talt_scaf_stop\tparent_start\tparent_stop\talt_start_tail\talt_stop_tail\n";



#---NUCMER TO CHARACTERIZE ALIGNMENTS---
my $job=0;
foreach my $htig(sort(keys(%output))){
    my $htig_out = $htig;
    $htig_out =~ s/\|.+//;       # remove |quiver|arrow|etc tags
    my $h = $output{$htig};
    if ($$h{top}){
        $available_threads->down(1);
        threads->create(\&nuc_align, $htig, $htig_out, $job);
        $job++;
    }
}

# wait on nucmer jobs
sleep 1 while threads->list();


#---WRITE THE PLACEMENT FILE---
for my $htig(sort(keys(%output))){
    my $h = $output{$htig};
    (($$h{hlen}) && ($$h{hlen} > (($coverage / 100) * $$h{len}))) ?                     # write the placement if > alignment cutoff to the best hit primary contig
        (print $OUT "Haplotigs\tPrimary_Assembly\t$htig\tSCAFFOLD\t$$h{top}\t",         # contig/assembly names
        ($$h{rlen} < 0 ? "-" : "+"),                                                    # check orientation (ncbi's 'mixed' tag not implemented)
        "\t$$h{hs}\t$$h{he}\t$$h{rs}\t$$h{re}\t",                                       # start/stop positions
        ($$h{hs} == 1 ? 0 : ($$h{hs}-1)), "\t",                                         # check start tail
        ($$h{he} == $$h{len} ? 0 : ($$h{len} - $$h{he})), "\n")                         # check end tail
    : print $OUT "Haplotigs\tPrimary_Assembly\t$htig", ("\tna" x 9), "\n";          # no placements if no or poor alignment
}

# DONE.
exit(0);


sub nuc_align {
    my $htig = shift;
    my $htig_out = shift;
    my $job = shift;
    
    my $h = $output{$htig};
    
    msg("processing haplotig $htig -> $$h{top}");
    # run nucmer etc if needed
    if (!(-s "$TMP/$htig_out.coords")){
        qruncmd("samtools faidx $haplotigs \"$htig\" > $TMP/$job.htig.fa");
        qruncmd("samtools faidx $primary \"$$h{top}\" > $TMP/$job.ref.fa");
        qruncmd("nucmer $TMP/$job.htig.fa $TMP/$job.ref.fa -p $TMP/$job.tmp");
        qruncmd("delta-filter -g $TMP/$job.tmp.delta > $TMP/$job.tmp.g.delta");
        qruncmd("show-coords -rT $TMP/$job.tmp.g.delta > $TMP/$htig_out.coords.tmp");
        qruncmd("mv $TMP/$htig_out.coords.tmp $TMP/$htig_out.coords");
    }

    # parse the coords file
    open my $COR, "$TMP/$htig_out.coords" or err("failed to open $TMP/$htig_out.coords");
    $$h{hs} = 999999999;
    $$h{rs} = 999999999;
    $$h{he} = 0;
    $$h{re} = 0;
    while(<$COR>){
        next if ($_ !~ /^\d/);
        my @l = split(/\s+/, $_);
        $$h{rlen} += $l[3] - $l[2];
        $$h{hlen} += $l[1] - $l[0];
        $$h{rs} < min($l[2],$l[3]) or $$h{rs} = min($l[2],$l[3]);
        $$h{re} > max($l[2],$l[3]) or $$h{re} = max($l[2],$l[3]);
        $$h{hs} < $l[0] or $$h{hs} = $l[0];
        $$h{he} > $l[1] or $$h{he} = $l[1];
    }
    close $COR;
    
    # cleanup and reset
    (unlink $_) for ("$TMP/$job.htig.fa", "$TMP/$job.ref.fa", "$TMP/$job.tmp.g.delta", "$TMP/$job.tmp.delta");
    
    $available_threads->up(1);
    threads->detach();
}





