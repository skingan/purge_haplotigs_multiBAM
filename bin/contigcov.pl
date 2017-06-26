#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# pipe script for analysing bedtools genomeCov output

my $incov;
my $outcsv;
my $lowc;
my $midc;
my $highc;
my $junk = 80;
my $suspect = 80;

my $usage = "
USAGE:
purge_haplotigs.pl contigcov  -i aligned.bam.genecov  -o coverage_stats.csv  -l 30  -m 80  -h 145  [ -j $junk  -s $suspect ]

REQUIRED:
-i      The bedtools genomecov output that was produced from 'purge_haplotigs readhist'
-o      Choose an output file name (csv format)
-l      The read depth low cutoff (use the histogram to eyeball these cutoffs)
-h      The read depth high cutoff
-m      The low point between the haploid and diploid peaks

OPTIONAL:
-j      Auto-assign contig as \"j\" (junk) if this percentage or greater of the contig is 
        low/high coverage (default = $junk, > 100 = don't junk anything)
-s      Auto-assign contig as \"s\" (suspected haplotig) if this percentage or less of the
        contig is diploid level of coverage (default = $suspect)

";


# globally declared for lazy subroutines
my $base_count;
my $n_low=0;
my $n_hap=0;
my $n_dip=0;
my $n_high=0;
my $low_p;
my $hap_p;
my $dip_p;
my $high_p;
my $b_hap_dip;
my $b_low_high;
my $cont_cur;
my $cont; 
my $cov;
my $bases;
my $all;
my $frac;

GetOptions(
    "in=s" => \$incov,
    "out=s" => \$outcsv,
    "low=i" => \$lowc,
    "high=i" => \$highc,
    "mid=i" => \$midc,
    "auto=i" => \$junk,
    "suspect=i" => \$suspect,
    "junk=i" => \$junk,
    "suspect=i" => \$suspect
) or die $usage;

if ( !($incov) || !($outcsv) || !($lowc) || !($highc) || !($midc) ){
    die $usage;
}

my $IN;
my $OUT;


open($IN, $incov) or die "failed to open $incov for reading\n";
open($OUT, ">$outcsv") or die "failed to open $outcsv for writing\n";

print $OUT "#contig,contig_reassign,bases_hap_dip,bases_low_high,bases_all,perc_low_coverage,perc_hap_coverage,perc_dip_coverage,perc_high_coverage\n";

while(<$IN>){
    ($cont, $cov, $bases, $all, $frac) = split(/\s+/, $_);
    # remove vert bar stuff like |quiver and |arrow if present
    $cont =~ s/\|.+//;
    # skip contigs with nothing
    next if ($cov == 0 && $frac == 1);
    # ignore genome histogram
    next if ($cont eq "genome");
    # initial iteration check
    $cont_cur = $cont if (!($cont_cur)); 
    if ($cont ne $cont_cur){
        print_contig();
                
        # reset for new contig
        $cont_cur = $cont;
        $n_low = 0; 
        $n_hap = 0; 
        $n_dip = 0; 
        $n_high = 0;
        
        # score
        collect_counts();
    } else {
        collect_counts();
    }
}
# print the final contig
print_contig();

close($IN);
close($OUT);

print STDERR "
Analysis finished successfully, contig coverage stats saved to $outcsv. You can modify this before
feeding it into the next step if you wish (i.e. if you wanted to add or remove 's' or 'j' flags 
for certain contigs).

";

exit(0);


### SUBROUTINES ###

sub collect_counts {
    $n_low += $bases if ($cov < $lowc);
    $n_hap += $bases if ($cov >= $lowc && $cov <= $midc);
    $n_dip += $bases if ($cov > $midc && $cov <= $highc);
    $n_high += $bases if ($cov > $highc);
}

sub print_contig {
    $base_count = $n_low + $n_hap + $n_dip + $n_high;
    
    $low_p = sprintf("%.3f", (100 * ($n_low / $base_count)));
    $hap_p = sprintf("%.3f", (100 * ($n_hap / $base_count)));
    $dip_p = sprintf("%.3f", (100 * ($n_dip / $base_count)));
    $high_p = sprintf("%.3f", (100 * ($n_high / $base_count)));
    
    $b_hap_dip = $n_hap + $n_dip;
    $b_low_high = $n_low + $n_high;
    
    my $assign = "";
    
    if ($junk){
        if ( ($low_p + $high_p) >= $junk ){
            $assign = "j";
        }
    }
    if ($suspect){
        if ($dip_p <= $suspect){
            $assign = "s" if (!($assign));
        }
    }
    
    print $OUT "$cont_cur,$assign,$b_hap_dip,$b_low_high,$base_count,$low_p,$hap_p,$dip_p,$high_p\n";
}


