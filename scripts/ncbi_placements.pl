#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PipeUtils;
use List::Util qw(min max);


# Adaptation of nucmer2ncbiPlacement.py from https://github.com/skingan/NCBI_DiploidAssembly/blob/master/nucmer2ncbiPlacement.py
# The adaptation is necessary if the assembly has been curated with purge_haplotigs as the FALCON IDs won't necessarily match up.


# PIPELINE:
    # blastn haplotigs against primary contigs and get best hit for each
    # nucmer, delta-filter, show-coords each haplotig against it's blastn hit
    # get the alignment information from the coords file and write to the output placements file


my $primary;
my $haplotigs;
my $out = "ncbi_placements.tsv";
my $threads = 4;
my $TMP = "temp_ncbi_placements";


my %output;     # $output{$htig}{len} = h-tig length
                #               {top} = top hit
                #               {rlen} = ref alignments running total (negative = rev alignments)
                #               {hlen} = htig alignments running total
                #               {hs} = htig start
                #               {he} = htig end
                #               {rs} = ref start
                #               {re} = ref end


my $usage = "
USAGE:
ncbi_placements.pl  -p primary_contigs.fasta  -h haplotigs.fasta  [ -o outfile  -t threads]

NOTE: designed for nucmer 3.9 or higher

REQUIRED:
-p      primary contigs fasta file
-h      haplotigs fasta file

OPTIONAL:
-o      out file name (DEFAULT = $out)
-t      threads for blastn and nucmer steps (DEFAULT = $threads)
";


# check dependencies
check_programs("blastn", "samtools", "parallel", "nucmer", "delta-filter", "show-coords");


# parse args
GetOptions (
    "primary=s" => \$primary,
    "haplotigs=s" => \$haplotigs,
    "out=s" => \$out,
    "threads=i" => \$threads
) or die $usage;


# check args and files
($primary) && ($haplotigs) || err("missing command line argument(s)\n$usage");

check_files($primary, $haplotigs);


# make the temp dir
mkdir $TMP if (!(-d $TMP));


# get fasta indexes ready and get the haplotig IDs and lengths
((-s "$_.fai") || runcmd({ command => "samtools faidx $_ 2> $TMP/samtools_faidx.stderr", logfile => "$TMP/samtools_faidx.stderr" })) for ($haplotigs, $primary);

# read in the contig lengths
open my $HFI, "<", "$haplotigs.fai" or err("failed to open file $haplotigs.fai");
while(<$HFI>){
    my @l = split(/\s+/, $_);
    $output{$l[0]}{len} = $l[1];
}
close $HFI;


# run blastn to get the best hit for each haplotig
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


# for each haplotig: run nucmer against top hit, delta-filter, and show-coords; parse the coords file to generate ncbi placement outputs
foreach my $htig(sort(keys(%output))){
    my $h = $output{$htig};
    if ($$h{top}){
        msg("processing haplotig $htig -> $$h{top}");
        # run nucmer etc if needed
        if (!(-s "$TMP/$htig.coords")){
            qruncmd("samtools faidx $haplotigs \"$htig\" > $TMP/htig.fa");
            qruncmd("samtools faidx $primary \"$$h{top}\" > $TMP/ref.fa");
            qruncmd("nucmer -t $threads $TMP/htig.fa $TMP/ref.fa -p $TMP/tmp");
            qruncmd("delta-filter -g $TMP/tmp.delta > $TMP/tmp.g.delta");
            qruncmd("show-coords -rT $TMP/tmp.g.delta > $TMP/$htig.coords");
        }
        
        # parse the coords file
        open my $COR, "$TMP/$htig.coords" or err("failed to open $TMP/$htig.coords");
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
    }
    
    (($$h{hlen}) && ($$h{hlen} > (0.5 * $$h{len}))) ?                               # write the placement if > 50% alignment to the best hit primary contig
        (print $OUT "Haplotigs\tPrimary_Assembly\t$htig\tSCAFFOLD\t$$h{top}\t",         # contig/assembly names
        ($$h{rlen} < 0 ? "-" : "+"),                                                    # check orientation (ncbi's 'mixed' tag not implemented)
        "\t$$h{hs}\t$$h{he}\t$$h{rs}\t$$h{re}\t",                                       # start/stop positions
        ($$h{hs} == 1 ? 0 : ($$h{hs}-1)), "\t",                                         # check start tail
        ($$h{he} == $$h{len} ? 0 : ($$h{len} - $$h{he})), "\n")                         # check end tail
    : print $OUT "Haplotigs\tPrimary_Assembly\t$htig", ("\tna" x 9), "\n";          # no placements if no or poor alignment
}


# DONE.
exit(0);




