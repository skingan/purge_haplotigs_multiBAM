#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use threads;
use Thread::Semaphore;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use PipeUtils;


#---SET UP LOGGING---
our $LOG;
my $TMP_DIR = "tmp_purge_haplotigs";
(-d $TMP_DIR) or mkdir $TMP_DIR;
open $LOG, ">", "$TMP_DIR/purge_haplotigs_purge.log" or die "failed to open log file for writing";



#---INPUT PARAMETERS---

my $genome;
my $coverage_stats;
my $bam_file;
my $threads = 4;
my $maxmatch_cutoff = 250;
my $bestmatch_cutoff = 70;
my $low_cutoff = 40;
my $out_prefix = "curated";
my $window_size = 9000;
my $step_size = 3000;
my $blastn_parameters = "-evalue 0.00001 -max_target_seqs 50 -max_hsps 1000 -word_size 24 -culling_limit 10";
my $lastz_parameters = "--gfextend --chain --gapped --seed=match14 --allocate:traceback=200.0M";


#---HELP MESSAGE---

my $usage = "
USAGE:
purge_haplotigs  purge  -g genome.fasta  -c coverage_stats.csv  -b aligned.bam

REQUIRED:
-g / -genome        Genome assembly in fasta format. Needs to be indexed with samtools faidx.
-c / -coverage      Contig by contig coverage stats csv file from the previous step.
-b / -bam           Samtools-indexed bam file of aligned and sorted reads/subreads to the reference 
                    (same one used for generating the read-depth histogram).

OPTIONAL:
-t / -threads       Number of worker threads to use. DEFAULT = $threads.
-o / -outprefix     Prefix for the curated assembly. DEFAULT = \"$out_prefix\".
-a / -align_cov     Percent cutoff for identifying a contig as a haplotig. DEFAULT = $bestmatch_cutoff.
-m / -max_match     Percent cutoff for identifying repetitive contigs. DEFAULT = $maxmatch_cutoff.
-wind_len           Length of genome window (in bp) for the coverage track in the dotplots. DEFAULT = $window_size.
-wind_step          Step distance for genome windows for coverage track. DEFAULT = $step_size.

";


#---CHECK PROGRAMS---

if (!(check_programs("bedtools", "blastn", "lastz", "makeblastdb", "samtools"))){
    err("ONE OR MORE REQUIRED PROGRAMS IS MISSING");
} else {
    msg("All required dependencies OK");
}

my $gnuparallel = check_programs("parallel");
if (!($gnuparallel)){
    print_message("WARN: not using GNU parallel, the blastn search will take longer to run");
}



#---PARSE ARGUMENTS---

my $args = "@ARGV";
GetOptions (
    "genome=s" => \$genome,
    "coverage_stats=s" => \$coverage_stats,
    "bam=s" => \$bam_file,
    "threads=i" => \$threads,
    "outprefix=s" => \$out_prefix,
    "align_cov=i" => \$bestmatch_cutoff,
    "max_match=i" => \$maxmatch_cutoff, 
    "wind_len=i" => \$window_size,
    "wind_step=i" => \$step_size
) or die $usage;

($genome) && ($coverage_stats) && ($bam_file) or die $usage;



# check the files
if (!(check_files($genome, $coverage_stats, $bam_file))){
    msg("one or more files missing, exiting");
    die $usage;
}



#---GLOBAL VARIABLES---

# files etc
my $MINCE_DIR = "$TMP_DIR/minced";
my $MINCE_DONE = "$TMP_DIR/minced.done";
my $BED_COVERAGE_DONE = "$TMP_DIR/bed_coverage.done";
my $LASTZ_DIR = "$TMP_DIR/tmp_lastz";
my $ASSIGNED = "dotplots_reassigned_contigs";
my $UNASSIGNED = "dotplots_unassigned_contigs";

my $tmp_asm = "$TMP_DIR/assembly.fasta";
my $suspects_fasta = "$TMP_DIR/suspects.fasta";
my $blastdb = "$TMP_DIR/blastdb/assembly";
my $blastgz = "$TMP_DIR/suspects.blastn.gz";
my $blast_summary = "$TMP_DIR/blastn_summary.tsv";
my $suspect_reassign = "$TMP_DIR/suspect_reassignments.tsv";
my $lastz_log = "$TMP_DIR/lastz_analysis.stderr";

# some clean file names
my $bam_file_name = $bam_file;
my $genome_file_name = $genome;
$bam_file_name =~ s/.*\///;
$genome_file_name =~ s/.*\///;

# reassignment step
my %suspects;   # suspects flagged from blastn search
my %junk;       # junk flagged from coverage analysis

my %contigs;    # $contigs{contig}{LEN} = 1 000 000     (length of contig)
                #                 {1} = contig          (first best blastn hit)
                #                 {2} = contig          (second best blastn hit)
                #                 {BM} = 95.0           (best match coverage percentage)
                #                 {MM} = 400.0          (max match coverage percentage)
                #                 {ASSIGN} = r/h/n/u    (contig reassignment, r=repeat, h=haplotig, n=no reassign, u=unknown)
                #                 {REASSIGNED} = 1/0    (if it's been reassigned)
                #                 {SUFFIX} = " <<HTIG"  (for reassigned ctgs, will be the fasta description)

# contig association step:
my %primaries;
my %refs;
my @current_path;
my @current_path_rename;
my $current_depth;



#---OUTPUT FILES---

my $contig_paths = "$out_prefix.contig_associations.log";
my $out_fasta = "$out_prefix.fasta";
my $out_haplotigs = "$out_prefix.haplotigs.fasta";
my $out_artefacts = "$out_prefix.artefacts.fasta";
my $out_reassignments = "$out_prefix.reassignments.tsv";



#---THREADS---

my $available_threads = Thread::Semaphore->new($threads);
my $writing_to_out = Thread::Semaphore->new(1);



#---FILEHANDLES---
my $FAI;
my $CSV;
my $GEN;
my $BLS;
my $TSV;
my $MAS;
my $MCV;
my $PTH;
my $CUH;
my $CUT;



#---OPEN DIRS ETC---

foreach my $dir ($TMP_DIR, $MINCE_DIR, $LASTZ_DIR, $ASSIGNED, $UNASSIGNED){
    if (!(-d $dir)){
        mkdir $dir or err("failed to create directory $dir");
    }
}


foreach my $file ($out_artefacts, $out_fasta, $out_haplotigs, $out_reassignments, $contig_paths){
    if (-s $file){
        unlink $file or err("failed to clean up previous run output file: $file");
    }
}



#---PRINT PARAMETERS---

msg("
Genome fasta:           $genome
Coverage csv:           $coverage_stats
Bam file:               $bam_file
Threads:                $threads
Cutoff, alignment:      $bestmatch_cutoff %
Cutoff, repeat:         $maxmatch_cutoff %
Cutoff, suspect:        $low_cutoff %
Out prefix:             $out_prefix
Coverage window len:    $window_size bp
Window step dist:       $step_size bp
Blastn parameters:      $blastn_parameters
Lastz paramters:        $lastz_parameters

Running using command:
purge_haplotigs purge $args\n
");



#---PIPELINE BEGIN---


msg("\n\n###\n\nBEGINNING PIPELINE!\n\n###\n");

# read in fasta.fai
read_fasta_fai();

# read in coverage stats
read_cov_stats();

# mince genome, this will make later steps run much faster
mince_genome();

# make bed windows, reads per window, for each contig. Used later for juxtaposing read-depth coverage with dotplots
get_window_coverage();

# make the blastdb and suspects.fasta, and run the blastn search
run_blastn();

# convert blast output to a simple list for quick iteration
hit_summary();


#---ITERATIVE STEP---

my $convergence = 0;
my $passes = 1;


while(!($convergence)){
    
    msg("\n\n###\n\nRUNNING PURGING PASS $passes\n\n###\n");
        
    # read through blastn hit summary and get top 2 matches for each suspect contig
    get_contig_hits();

    # run mummer steps
    run_lastz_analysis();
    
    # remove conflict reassignments
    check_assignments();
    
    # add to reassignments list
    add_reassignments();
}


#---GENERATE OUTPUT---


msg("\n\n###\n\nGENERATING OUTPUT\n\n###\n");

# get the reassignment paths
get_reassignment_paths();

# write the table and new assembly
write_assembly();


msg("\n\n###\n\nPURGE HAPLOTIGS HAS COMPLETED SUCCESSFULLY!\n\n###\n");



#---SUBROUTINES---


sub read_fasta_fai {
    (-s "$genome.fai") or runcmd({ command => "samtools faidx $genome 2> $TMP_DIR/samtools.faidx.stderr", logfile => "$TMP_DIR/samtools.faidx.stderr" });
    msg("reading in genome index file: $genome.fai");
    
    open $FAI, "<",  "$genome.fai" or err ("failed to open  $genome.fai for reading");
    
    while(my $l = <$FAI>){
        my @line = split(/\s+/, $l);
        $line[0] =~ s/\|.*$//;
        if (!($line[0]) || !($line[1])){
            err("bad entry in $genome.fai index file? line:\n$l");
        } else {
            $contigs{$line[0]}{LEN} = $line[1];
        }
    }
    
    close $FAI;
}



sub read_cov_stats {
    msg("reading in contig coverage stats file: $coverage_stats");
    
    open $CSV, "<", $coverage_stats or err("failed to open $coverage_stats for reading");
    
    while(my $l = <$CSV>){
        next if ($l =~ /^#/);
        
        my @line = split(/,/, $l);
        if (($line[0]) && !($contigs{$line[0]})){
            err("no contig \"$line[0]\" in genome index file, csv line:\n$l");
        } elsif ($line[1] =~ /^s$/i){
            $suspects{$line[0]} = 1;
        } elsif ($line[1] =~ /^j$/i){
            $junk{$line[0]} = 1;
        }
    }
}



sub mince_genome {
    if ( (-e $MINCE_DONE) && (-d $MINCE_DIR) ){
        msg("skipping already completed step \"mince genome\". If you wanted to re-run this step please delete \"$MINCE_DONE\" and rerun this pipeline");
    } else {
        msg("\"mincing\" genome");
        open $GEN, "<", $genome or err("failed to open $genome for reading");
        my $OUT;
        
        while (my $l = <$GEN>) {
            if ($l !~ /^>/) {
                print $OUT $l;
            } else {
                my $id;
                if ($l =~ /^>([a-zA-Z0-9_-]+)[\s\|]/){ # cut off |quiver|arrow|etc
                    $id = $1;
                    close $OUT if ($OUT);
                    open $OUT, ">", "$MINCE_DIR/$id.fasta" or err("failed to open \"$MINCE_DIR/$id.fasta\" for writing");
                    print $OUT $l;
                } else {
                    err("failed to get id from fasta file");
                }
            }
        }
        
        close $OUT;
        close $GEN;
        
        qruncmd("touch $MINCE_DONE");
    }
}



sub get_window_coverage {
    if (-e $BED_COVERAGE_DONE){
        msg("Skipping already completed step 'get window coverage'. If you want to rerun this step delete \"$BED_COVERAGE_DONE\" and rerun this pipeline");
    } else {
        msg("Getting windowed read depth coverage for each contig");
        
        (-s "$bam_file.bai") or runcmd({ command => "samtools index $bam_file 2> $TMP_DIR/samtools.index.sdterr", logfile => "$TMP_DIR/samtools.index.sdterr" });
        
        # make the bed windows
        runcmd({ command => "bedtools makewindows -g $genome.fai -w $window_size -s $step_size > $TMP_DIR/$genome_file_name.windows.bed.tmp 2> $TMP_DIR/bedtools.mkwind.stderr  &&  mv $TMP_DIR/$genome_file_name.windows.bed.tmp $TMP_DIR/$genome_file_name.windows.bed", logfile => "$TMP_DIR/bedtools.mkwind.stderr" });
        
        # get the reads per window
        runcmd({ command => "bedtools multicov -bams $bam_file -bed $TMP_DIR/$genome_file_name.windows.bed > $TMP_DIR/$genome_file_name.mcov.tmp 2> $TMP_DIR/bedtools.mcov.stderr  &&  mv $TMP_DIR/$genome_file_name.mcov.tmp $TMP_DIR/$genome_file_name.mcov", logfile => "$TMP_DIR/bedtools.mcov.stderr" });
        
        # generate histogram table files for each contig, dump them in the minced directory
        open $MCV, "$TMP_DIR/$genome_file_name.mcov" or err("Failed to open $TMP_DIR/$genome_file_name.mcov for reading");
        
        my @wind_cov_out;
        my $cur_ctg;
        
        while(<$MCV>){
            $_ =~ s/\n//;
            if ($_){
                my @l = split(/\s/, $_);
                $l[0] =~ s/\|.+//; # cut off |quiver|arrow|etc
                if (!($cur_ctg)){
                    $cur_ctg = $l[0];
                } elsif ($cur_ctg ne $l[0]){
                    # dump the output when we hit a new contig
                    write_cov_out($cur_ctg, \@wind_cov_out);
                    # reset
                    $cur_ctg = $l[0];
                    @wind_cov_out = ();
                }
                my $p = ($l[1] + $l[2]) / 2;
                my $d = $l[3] / ($l[2] - $l[1]);
                push @wind_cov_out, "$p\t$d\n";
            }
        }
        # dump the last contig
        write_cov_out($cur_ctg, \@wind_cov_out);
        
        sub write_cov_out {
            my $ctg = $_[0];
            my @wind_cov = @{$_[1]};
            open my $TMC, ">$MINCE_DIR/$ctg.cov" or err("Failed to open $MINCE_DIR/$ctg.cov for writing");
            print $TMC @wind_cov;
            close($TMC);
        }
        
        # touch the done file
        qruncmd("touch $BED_COVERAGE_DONE");
    }
}



sub run_blastn {
    if (-s $blastgz){
        msg("skipping already completed blastn step. To rerun this step please delete $blastgz and rerun this pipeline");
    } else {
        
        msg("preparing blastdb");
        
        if (-s $tmp_asm){
            unlink $tmp_asm or err("failed to clean up temp file $tmp_asm");
        }
        if (-s $suspects_fasta){
            unlink $suspects_fasta or err("failed to clean up temp file $suspects_fasta");
        }
        
        foreach my $contig (sort(keys(%contigs))){
            if (!($junk{$contig})){
                qruncmd("cat $MINCE_DIR/$contig.fasta >> $tmp_asm");
            }
            if ($suspects{$contig}){
                qruncmd("cat $MINCE_DIR/$contig.fasta >> $suspects_fasta");
            }
        }
        
        runcmd({ command => "makeblastdb -in $tmp_asm -dbtype nucl -out $blastdb 2>&1 1> $TMP_DIR/makeblastdb.stderr", logfile => "$TMP_DIR/makeblastdb.stderr" });
        
        msg("running blastn analysis");
        
        if ($gnuparallel){
            runcmd({ command => "cat $suspects_fasta | parallel --no-notice -j $threads --block 50k --recstart '>' --pipe blastn -db $blastdb -outfmt 6 $blastn_parameters -query - 2> $TMP_DIR/gnuparallel.stderr | awk ' \$1 != \$2 && \$4 > 500 { print } ' | gzip - > $blastgz.tmp", logfile => "$TMP_DIR/gnuparallel.stderr" });
        } else {
            runcmd({ command => "blastn -query $suspects_fasta -db $blastdb -outfmt 6 -num_threads $threads $blastn_parameters 2> $TMP_DIR/blastn.stderr |  awk ' \$1 != \$2 && \$4 > 500 { print } ' | gzip - > $blastgz.tmp.gz", logfile => "$TMP_DIR/blastn.stderr" });
        }
        qruncmd("mv $blastgz.tmp $blastgz");
    }
}



sub hit_summary {
    msg("preparing blastn hit summary");
    
    if (-s $blast_summary){
        unlink $blast_summary or err("failed to clean up temp file $blast_summary");
    }
    
    open $BLS, "gunzip -c $blastgz |" or err("failed to open gunzip pipe from $blastgz");
    
    open $TSV, ">", "$blast_summary.tmp" or err("failed to open $blast_summary.tmp for writing");
    
    my %uniq;
    
    while(my $l = <$BLS>){
        $l =~ s/\|\S+\s/\t/g;
        my @line = split(/\s+/, $l);
        next if ($line[0] eq $line[1]);
        if (!($uniq{$line[0].$line[1]})){
            print $TSV "$line[0]\t$line[1]\n";
            $uniq{$line[0].$line[1]} = 1;
        }
    }
    
    undef %uniq;
    close $BLS;
    close $TSV;
    
    qruncmd("mv $blast_summary.tmp $blast_summary");
}


#---ITERATIVE STEP---


sub get_contig_hits {
    msg("getting contig hits from blastn output");
    
    # remove reference contigs if they themselves have been reassigned 
    foreach my $ctg (sort(keys(%contigs))){
        next if ($contigs{$ctg}{REASSIGNED});
        if ($contigs{$ctg}{1}){
            if ($contigs{$contigs{$ctg}{1}}{REASSIGNED}){
                undef $contigs{$ctg}{1};
                undef $contigs{$ctg}{2};
                $contigs{$ctg}{ASSIGN} = 0;
            }
        }
        if ($contigs{$ctg}{2}){
            if ($contigs{$contigs{$ctg}{2}}{REASSIGNED}){
                undef $contigs{$ctg}{2};
                $contigs{$ctg}{ASSIGN} = 0;
            }
        }
    }
    
    # add best non-reassigned reference contigs for each contig 
    open $TSV, $blast_summary or err("failed to open $blast_summary for reading");
    
    while(my $l = <$TSV>){
        my @line = split(/\s+/, $l);
        if (!($contigs{$line[0]}{REASSIGNED}) && !($contigs{$line[1]}{REASSIGNED})){
            if (!($contigs{$line[0]}{1})){
                $contigs{$line[0]}{1} = $line[1];
                $contigs{$line[0]}{ASSIGN} = 0;
            } elsif (!($contigs{$line[0]}{2}) || ($line[1] ne $contigs{$line[0]}{1})){
                $contigs{$line[0]}{2} = $line[1];
                $contigs{$line[0]}{ASSIGN} = 0;
            }
        }
    }
    
    close $TSV;
}



sub run_lastz_analysis {
    msg("Running lastz analysis on blastn hits");
    
    my $jobnumber = 0;
    
    if (-s $suspect_reassign){
        qruncmd("mv $suspect_reassign $suspect_reassign.$passes");
    }
    
    if (-s $lastz_log){
        unlink $lastz_log or err("failed to clean up temp file $lastz_log");
    }
    
    CTG: foreach my $contig (sort(keys(%contigs))){
        next CTG if ($contigs{$contig}{REASSIGNED});
        next CTG if ($contigs{$contig}{ASSIGN});
        
        # skip if no hits
        if (!($contigs{$contig}{1})){
            $contigs{$contig}{ASSIGN} = "n";
            if (-s "$UNASSIGNED/$contig.png"){
                unlink "$UNASSIGNED/$contig.png";
            }
            next CTG;
        }
        
        # run the lastz job
        $available_threads->down(1);
        threads->create(\&lastz_job, $jobnumber, $contig);
        $jobnumber++;
    }
    
    # wait on remaining jobs
    sleep 1 while threads->list();
}



# This could do with some refactoring
sub lastz_job {
    my $tmp_log;
    my $cmd;
    my $lastz_fail=0;

    my $job = $_[0];
    my $contig = $_[1];
    my $ref1 = $contigs{$contig}{1};
    my $ref2 = $contigs{$contig}{2} if ($contigs{$contig}{2});

    # run lastz against ref1
    $cmd = "lastz $lastz_parameters --format=general --rdotplot=$LASTZ_DIR/$job.1.rdotplot $MINCE_DIR/$contig.fasta $MINCE_DIR/$ref1.fasta > $LASTZ_DIR/$job.coords 2> $LASTZ_DIR/$job.lastz.stderr\n";
    $tmp_log .= "RUNNING: $cmd";
    runcmd({ command => $cmd, logfile => "$LASTZ_DIR/$job.lastz.stderr", silent => 1 });
    
    if (!(-s "$LASTZ_DIR/$job.1.rdotplot")){
        $ref1 = 0;
    }

    # run lastz against ref2
    if ($ref2){
        $cmd = "lastz $lastz_parameters --format=general --rdotplot=$LASTZ_DIR/$job.2.rdotplot $MINCE_DIR/$contig.fasta $MINCE_DIR/$ref2.fasta >> $LASTZ_DIR/$job.coords 2>> $LASTZ_DIR/$job.lastz.stderr\n";
        $tmp_log .= "RUNNING: $cmd";
        runcmd({ command => $cmd, logfile => "$LASTZ_DIR/$job.lastz.stderr", silent => 1 });
        if (!(-s "$LASTZ_DIR/$job.2.rdotplot")){
            $ref2 = 0;
        }
    }

    # capture lastz STDERR message
    $tmp_log .= `cat $LASTZ_DIR/$job.lastz.stderr`;

    # sort the lastz alignments
    qruncmd("grep -v -P \"^#\" $LASTZ_DIR/$job.coords | sort -k5,5n -k6,6n > $LASTZ_DIR/$job.s.coords");

    # get 'maxmatch' and 'bestmatch' coverages from the sorted general output of lastz
    open my $BMC, "$LASTZ_DIR/$job.s.coords" or err("failed to open grep/sort pipe from $LASTZ_DIR/$job.s.coords for reading");

    my @p;
    my $bestmatch=0;
    my $maxmatch=0;
    while(<$BMC>){
        next if ($_ =~ /^#/);
        $_ =~ s/^\s+//;
        my @l = split(/\s+/, $_);
        $maxmatch+=($l[5]-$l[4]);
        if (@p){
            if ($l[4] > $p[1]){
                $bestmatch+=($p[1]-$p[0]);
                @p=($l[4], $l[5]);
            } elsif ($p[1] < $l[5]) {
                $p[1] = $l[5];
            } 
        } else {
            @p=($l[4], $l[5]);
        }
    }
    $bestmatch+=($p[1]-$p[0]) if (@p);
    close $BMC;
    $maxmatch = sprintf "%.2f", ($maxmatch/$contigs{$contig}{LEN}) * 100;
    $bestmatch = sprintf "%.2f", ($bestmatch/$contigs{$contig}{LEN}) * 100;

    $tmp_log .= "BESTMATCH coverage = $bestmatch\n";
    $tmp_log .= "MAXMATCH coverage = $maxmatch\n";

    my $assignment;

    # guess the assignment
    if ($bestmatch >= $bestmatch_cutoff){
        if ($maxmatch >= $maxmatch_cutoff){
            $assignment = "r";
        } else {
            $assignment = "h";
        }
    } elsif ($bestmatch < $low_cutoff){
        $assignment = "n";
    } else {
        $assignment = "u";
    }

    # make the dotplot in unassigned dir
    if (-s "$UNASSIGNED/$contig.png"){
        unlink "$UNASSIGNED/$contig.png" or err("failed to clean up previous dotplot $UNASSIGNED/$contig.png");
    }

    # perform a check to see if either or both reference alignments come up empty
    undef $cmd;
    if ($assignment ne "n"){
        if (($ref2) && ($ref1)){
            $cmd = "$RealBin/../scripts/dot_plot_plus.Rscript $UNASSIGNED/$contig.png $contig $contigs{$contig}{LEN} $MINCE_DIR/$contig.cov $ref1 $LASTZ_DIR/$job.1.rdotplot $ref2 $LASTZ_DIR/$job.2.rdotplot 1> $LASTZ_DIR/$job.Rscript.stderr 2>&1\n";
        } elsif (($ref1) && !($ref2)) {
            $cmd = "$RealBin/../scripts/dot_plot_plus.Rscript $UNASSIGNED/$contig.png $contig $contigs{$contig}{LEN} $MINCE_DIR/$contig.cov $ref1 $LASTZ_DIR/$job.1.rdotplot 1> $LASTZ_DIR/$job.Rscript.stderr 2>&1\n";
        } elsif (($ref2) && !($ref1)){
            $cmd = "$RealBin/../scripts/dot_plot_plus.Rscript $UNASSIGNED/$contig.png $contig $contigs{$contig}{LEN} $MINCE_DIR/$contig.cov $ref2 $LASTZ_DIR/$job.2.rdotplot 1> $LASTZ_DIR/$job.Rscript.stderr 2>&1\n";
        }

        if ($cmd){
            $tmp_log .= $cmd;
            runcmd({ command => $cmd, logfile => "$LASTZ_DIR/$job.Rscript.stderr", silent => 1 });
        } else {
            err("Contig $contig returned alignment score but no rdotplot files");
        }
    }

    # print the lastz log
    print_lastz_log(\$tmp_log);

    # print the reassignment
    $writing_to_out->down(1);
    open my $TSV, ">>", $suspect_reassign or err("failed to open $suspect_reassign for appending");
    if ($contigs{$contig}{2}){
        print $TSV "$contig\t$contigs{$contig}{1}\t$contigs{$contig}{2}\t$bestmatch\t$maxmatch\t$assignment\n";
    } else {
        print $TSV "$contig\t$contigs{$contig}{1}\t\t$bestmatch\t$maxmatch\t$assignment\n";
    }
    close $TSV;
    $writing_to_out->up(1);

    # clean up 
    foreach my $file ("$LASTZ_DIR/$job.lastz.stderr", "$LASTZ_DIR/$job.coords", "$LASTZ_DIR/$job.s.coords", "$LASTZ_DIR/$job.1.rdotplot", "$LASTZ_DIR/$job.2.rdotplot", "$LASTZ_DIR/$job.Rscript.stderr"){
        if (-e $file){
            unlink $file or err("failed to clean up temp file $file");
        }
    }

    # exit
    $available_threads->up(1);
    threads->detach();


    sub print_lastz_log {
        $writing_to_out->down(1);
        open my $LLOG, ">>", $lastz_log or err("failed to open $lastz_log for appended writing");
        print $LLOG ${$_[0]};
        close $LLOG;
        $writing_to_out->up(1);
    }
}



sub check_assignments {
    msg("Checking contig assignments for conflicts");
    
    # read in all the mummer assignments
    if (-s $suspect_reassign){
        open $MAS, $suspect_reassign or err("failed to open $suspect_reassign for reading");
        
        while(my $l = <$MAS>){
            my @line = split(/\t/, $l);
            $line[5] =~ s/\s//g;
            $contigs{$line[0]}{ASSIGN} = $line[5];
            $contigs{$line[0]}{BM} = $line[3];
            $contigs{$line[0]}{MM} = $line[4];
        }
        
        close $MAS;
    }
    
    # check all assignments for conflicts
    foreach my $ctg (sort(keys(%contigs))){
        next if ($contigs{$ctg}{ASSIGN} !~ /[rh]/i);
        next if ($contigs{$ctg}{REASSIGNED});
        if ($contigs{$contigs{$ctg}{1}}{REASSIGNED}){
            err("ref seq was already reassigned, this should not have happened");
        }
        
        my $r_ctg = $contigs{$ctg}{1};
        if ($contigs{$r_ctg}{ASSIGN} =~ /[rh]/i){
            msg("conflict: $ctg and it's match $r_ctg both flagged for reassignment");
            
            # if both haplotigs, keep the longest
            if ( ($contigs{$ctg}{ASSIGN} =~ /h/i) && ($contigs{$r_ctg}{ASSIGN} =~ /h/i) ){
                if ($contigs{$ctg}{LEN} > $contigs{$r_ctg}{LEN}){
                    $contigs{$ctg}{ASSIGN} = 0;
                    msg("\tkeeping longer contig $ctg");
                } else {
                    $contigs{$r_ctg}{ASSIGN} = 0;
                    msg("\tkeeping longer contig $r_ctg");
                }
            }
            
            # if one is repetitive, keep haplotig
            elsif ( ($contigs{$ctg}{ASSIGN} . $contigs{$r_ctg}{ASSIGN}) =~ /rh|hr/i){
                if ($contigs{$ctg}{ASSIGN} =~ /h/i){
                    msg("\tkeeping haplotig $ctg, removing repeat $r_ctg");
                    $contigs{$ctg}{ASSIGN} = 0;
                } else {
                    msg("\tkeeping haplotig $r_ctg, removing repeat $ctg");
                    $contigs{$r_ctg}{ASSIGN} = 0;
                }
            }
            
            # if both repeats, keep the longest
            elsif ( ($contigs{$ctg}{ASSIGN} =~ /r/i) && ($contigs{$r_ctg}{ASSIGN} =~ /r/i) ){
                if ($contigs{$ctg}{LEN} > $contigs{$r_ctg}{LEN}){
                    msg("\tkeeping longer repeat contig $ctg, removing repeat $r_ctg");
                    $contigs{$ctg}{ASSIGN} = 0;
                } else {
                    msg("\tkeeping longer repeat contig $r_ctg, removing repeat $ctg");
                    $contigs{$r_ctg}{ASSIGN} = 0;
                }
            }
            
            else {
                err("Unknown combination, this should not have occured");
            }
        }
    }
}



sub add_reassignments {
    msg("Logging reassignments and checking for convergence");
    
    my $convergence_check = 1;
    
    foreach my $ctg(sort(keys(%contigs))){
        next if ($contigs{$ctg}{REASSIGNED});
        
        if ($contigs{$ctg}{ASSIGN} =~ /[rh]/i){
            $contigs{$ctg}{REASSIGNED} = 1;
            $convergence_check = 0;
            qruncmd("mv $UNASSIGNED/$ctg.png $ASSIGNED/$ctg.png");
        }
    }
    
    # convergence check
    if ($convergence_check){
        $convergence = 1;
        msg("convergence reached!");
    } else {
        msg("convergence not reached, more passes needed");
    }
    
    $passes++;
}


#---END ITERATIVE STEP---


sub get_reassignment_paths {
    msg("writing contig association paths");
    
    # get contigs that are references, and contigs that are references AND not themselves reassigned
    foreach my $ctg (sort(keys(%contigs))){
        if ($contigs{$ctg}{REASSIGNED}){
            push @{$refs{$contigs{$ctg}{1}}}, $ctg;
            if ( !($primaries{$contigs{$ctg}{1}}) && !($contigs{$contigs{$ctg}{1}}{REASSIGNED}) ){
                $primaries{$contigs{$ctg}{1}} = 1;
            }
            
            # rename the ASSIGN key for graph output
            if ($contigs{$ctg}{ASSIGN} eq "h"){
                $contigs{$ctg}{ASSIGN} = "HAPLOTIG";
            } elsif ($contigs{$ctg}{ASSIGN} eq "r"){
                $contigs{$ctg}{ASSIGN} = "REPEAT";
            } else {
                err("unknown reassignment $contigs{$ctg}{ASSIGN}, this shouldn't have happened")
            }
        }
    }
    
    open $PTH, ">", $contig_paths or err("failed to open $contig_paths for writing");
    
    # use the primaries to recursively find the association path
    foreach my $ctg (sort(keys(%primaries))){
        @current_path = ("$ctg,PRIMARY");
        @current_path_rename = @current_path;
        $current_depth = 0;
        find_path($ctg);
        print $PTH "\n";
    }
}



sub find_path {
    my $ref_ctg = $_[0];
    
    if (!($refs{$ref_ctg})){
        print $PTH @current_path;
        print $PTH "\n";
        my $suffix;
        foreach(@current_path_rename){
            my @ctg = split(/,/,$_);
            if (!($ctg[1] =~ /PRIMARY/)){
                err(@current_path_rename) if (!($suffix));
                $contigs{$ctg[0]}{REASSIGNED} = "$ctg[1]$suffix";
                $suffix = "<--$ctg[0]_$ctg[1]$suffix";
            } else {
                $suffix = "<--$ctg[0]_$ctg[1]";
            }
        }
    } else {
        foreach my $ctg (@{$refs{$ref_ctg}}){
            push @current_path, " -> $ctg,$contigs{$ctg}{ASSIGN}";
            push @current_path_rename, "$ctg,$contigs{$ctg}{ASSIGN}";
            $current_depth += 1;
            
            find_path($ctg);
            
            pop @current_path;
            pop @current_path_rename;
            $current_depth -= 1;
            for (my $i=0; $i<=$current_depth; $i++){
                $current_path[$i] =~ s/./ /g;
            }
        }
    }
}



sub write_assembly {
    msg("writing the reassignment table and new assembly files");
    
    open $CUH, ">", $out_haplotigs or err("failed to open $out_haplotigs for writing");
    open $CUT, ">", $out_reassignments or err("failed to open $out_reassignments for writing");
    
    print $CUT "#reassigned_contig\ttop_hit_contig\tsecond_hit_contig\tbest_match_coverage\tmax_match_coverage\treassignment\n";
    foreach my $ctg (sort(keys(%contigs))){
        if (($contigs{$ctg}{REASSIGNED}) && !($junk{$ctg})){
            write_seq($CUH, $ctg);
            my $c2 = $contigs{$ctg}{2} || "-";
            print $CUT "$ctg\t$contigs{$ctg}{1}\t$c2\t$contigs{$ctg}{BM}\t$contigs{$ctg}{MM}\t$contigs{$ctg}{ASSIGN}\n";
        }
    }
    
    print $CUT "#contigs_kept\n";
    foreach my $ctg (sort(keys(%contigs))){
        if ( !($junk{$ctg}) && !($contigs{$ctg}{REASSIGNED}) ){
            my $c1 = $contigs{$ctg}{1} || "-";
            my $c2 = $contigs{$ctg}{2} || "-";
            my $b = $contigs{$ctg}{BM} || "-";
            my $m = $contigs{$ctg}{MM} || "-";
            print $CUT "$ctg\t$c1\t$c2\t$b\t$m\tKEEP\n";
            qruncmd("cat $MINCE_DIR/$ctg.fasta >> $out_fasta");
        }
    }
    
    print $CUT "#junk_contigs\n";
    foreach my $ctg (sort(keys(%contigs))){
        if ($junk{$ctg}){
            qruncmd("cat $MINCE_DIR/$ctg.fasta >> $out_artefacts");
            print $CUT "$ctg\t-\t-\t-\t-\tJUNK\n";
        }
    }
}



sub write_seq {
    my $fh = $_[0];
    my $ctg = $_[1];
    
    open my $FA, "$MINCE_DIR/$ctg.fasta" or err("failed to open $MINCE_DIR/$ctg.fasta for reading");
    
    while(my $l = <$FA>){
        if ($l =~ /^>/){
            print $fh ">$ctg $contigs{$ctg}{REASSIGNED}\n";
        } else {
            print $fh $l;
        }
    }
}




