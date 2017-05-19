#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use threads;
use Thread::Semaphore;
use FindBin qw($Bin);



#---INPUT PARAMETERS---

my $genome;
my $coverage_stats;
my $threads = 4;
my $maxmatch_cutoff = 250;
my $bestmatch_cutoff = 75;
my $low_cutoff = 40;
my $out_prefix = "curated";



#---HELP MESSAGE---

my $usage = "
USAGE:

purge_haplotigs.pl  -genome genome.fasta  -coverage coverage_stats.csv

REQUIRED:
-g / -genome        Genome assembly in fasta format. Needs to be indexed with samtools faidx.
-c / -coverage      Contig by contig coverage stats csv file from the previous step.

OPTIONAL:
-t / -threads       Number of worker threads to use. DEFAULT = $threads.
-o / -outprefix     Prefix for the curated assembly. DEFAULT = \"$out_prefix\".
-b / -best_match    Percent cutoff for identifying a contig as a haplotig. DEFAULT = $bestmatch_cutoff.
-m / -max_match     Percent cutoff for identifying repetitive contigs. DEFAULT = $maxmatch_cutoff.

";


#---CHECK PROGRAMS---

my $dependencies = 1;
$dependencies x= chkprog("blastn", "-h > /dev/null 2>&1");
$dependencies x= chkprog("makeblastdb", "-h > /dev/null 2>&1");
$dependencies x= chkprog("$Bin/../mummer/nucmer", "-h > /dev/null 2>&1");
$dependencies x= chkprog("$Bin/../mummer/delta-filter", "-h > /dev/null 2>&1");
$dependencies x= chkprog("$Bin/../mummer/show-coords", "-h > /dev/null 2>&1");
$dependencies x= chkprog("$Bin/../mummer/mummerplot", "-h > /dev/null 2>&1");

if (!($dependencies)){
    err("ONE OR MORE REQUIRED PROGRAMS IS MISSING");
} else {
    msg("ALL DEPENDENCIES OK");
}

my $gnuparallel = chkprog("parallel", "--version > /dev/null 2>&1");
if (!($gnuparallel)){
    print_message("WARN: not using gnuparallel, the blastn search will take longer to run");
}



#---PARSE ARGUMENTS---

GetOptions (
    "genome=s" => \$genome,
    "coverage_stats=s" => \$coverage_stats,
    "threads=i" => \$threads,
    "outprefix=s" => \$out_prefix,
    "best_match=i" => \$bestmatch_cutoff,
    "max_match=i" => \$maxmatch_cutoff
) or die($usage);

die $usage if (!($genome) || !($coverage_stats));

if (!(check_files($genome, "$genome.fai", $coverage_stats))){
    msg("one or more files missing, exiting");
    die $usage;
}



#---GLOBAL VARIABLES---

# files etc
my $TMP_DIR = "tmp_purge_haplotigs";
my $MINCE_DIR = "$TMP_DIR/minced";
my $MINCE_DONE = "$TMP_DIR/minced.done";
my $MUM_DIR = "$TMP_DIR/tmp_mummer";
my $ASSIGNED = "dotplots_reassigned_contigs";
my $UNASSIGNED = "dotplots_unassigned_contigs";

my $tmp_asm = "$TMP_DIR/assembly.fasta";
my $suspects_fasta = "$TMP_DIR/suspects.fasta";
my $blastdb = "$TMP_DIR/blastdb/assembly";
my $blastgz = "$TMP_DIR/suspects.blastn.gz";
my $blastn_done = "$TMP_DIR/blastn.done";
my $blast_summary = "$TMP_DIR/blastn_summary.tsv";
my $suspect_reassign = "$TMP_DIR/suspect_reassignments.tsv";
my $mummer_log = "$TMP_DIR/mummer_analysis.stderr";

# reassignment step
my %suspects;   #   suspects flagged from blastn search
my %junk;       #   junk flagged from coverage analysis

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

my $LOG;
my $FAI;
my $CSV;
my $GEN;
my $BLS;
my $TSV;
my $MAS;

my $PTH;
my $CUH;
my $CUT;



#---OPEN DIRS ETC---

foreach my $dir ($TMP_DIR, $MINCE_DIR, $MUM_DIR, $ASSIGNED, $UNASSIGNED){
    if (!(-d $dir)){
        mkdir $dir or err("failed to create directory $dir");
    }
}

open $LOG, ">", "$TMP_DIR/purge_haplotigs.log" or err("failed to open $TMP_DIR/purge_haplotigs.log for writing");

foreach my $file ($out_artefacts, $out_fasta, $out_haplotigs, $out_reassignments, $contig_paths){
    if (-s $file){
        unlink $file or err("failed to clean up previous run output file: $file");
    }
}




#---PIPELINE BEGIN---

msg("\n\n###\n\nBEGINNING PIPELINE!\n\n###\n");

# read in fasta.fai
read_fasta_fai();

# read in coverage stats
read_cov_stats();

# mince genome, this will make later steps run much faster
mince_genome();

# make the blastdb and suspects.fasta, and run the blastn search
run_blastn();

# summarise hits
hit_summary();


#---ITERATIVE STEP---

my $convergence = 0;
my $passes = 1;

while(!($convergence)){
    
    msg("\n\n###\n\nRUNNING PURGING PASS $passes\n\n###\n");
        
    # read through blastn hit summary and get top 2 matches for each suspect contig
    get_contig_hits();

    # run mummer steps
    run_mummer_analysis();
    
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
                if ($l =~ /^>([a-zA-Z0-9_-]+)[\s\|]/){
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

sub run_blastn {
    if ( (-e $blastn_done) && (-s $blastgz) ){
        msg("skipping already completed blastn step. To rerun this step please delete $blastn_done and rerun this pipeline");
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
        
        runcmd("makeblastdb -in $tmp_asm -dbtype nucl -out $blastdb 2>&1 1>", "$TMP_DIR/makeblastdb.stderr");
        
        msg("running blastn analysis");
        
        if ($gnuparallel){
            runcmd("cat $suspects_fasta | parallel --no-notice -j $threads --block 100k --recstart '>' --pipe blastn -db $blastdb -outfmt 6 -evalue 0.000000000001 -max_target_seqs 20 -max_hsps 1000 -word_size 28 -culling_limit 10 -query - 2> $TMP_DIR/gnuparallel.stderr | awk ' \$1 != \$2 && \$4 > 500 { print } ' | gzip - > $blastgz");
        } else {
            runcmd("blastn -query $suspects_fasta -db $blastdb -outfmt 6 -evalue 0.000000000001 -num_threads $threads -max_target_seqs 3 -max_hsps 1000 -word_size 28 -culling_limit 10 2> $TMP_DIR/blastn.stderr |  awk ' \$1 != \$2 && \$4 > 500 { print } ' | gzip - > $blastgz");
        }
        
        qruncmd("touch $blastn_done");
    }
}

sub hit_summary {
    msg("preparing blastn hit summary");
    
    if (-s $blast_summary){
        unlink $blast_summary or err("failed to clean up temp file $blast_summary");
    }
    
    open $BLS, "gunzip -c $blastgz |" or err("failed to open gunzip pipe from $blastgz");
    
    open $TSV, ">", $blast_summary or err("failed to open $blast_summary for writing");
    
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
}

#---ITERATIVE STEP---

sub get_contig_hits {
    msg("getting contig hits from blastn output");
    
    foreach my $ctg (sort(keys(%contigs))){
        next if ($contigs{$ctg}{REASSIGNED});
        if ($contigs{$ctg}{1}){
            if ($contigs{$contigs{$ctg}{1}}{REASSIGNED}){
                undef $contigs{$ctg}{1};
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
    
    open $TSV, $blast_summary or err("failed to open $blast_summary for reading");
    
    while(my $l = <$TSV>){
        my @line = split(/\s+/, $l);
        if (!($contigs{$line[0]}{REASSIGNED}) && !($contigs{$line[1]}{REASSIGNED})){
            if (!($contigs{$line[0]}{1})){
                $contigs{$line[0]}{1} = $line[1];
            } elsif (!($contigs{$line[0]}{2})){
                $contigs{$line[0]}{2} = $line[1];
            }
        }
    }
    
    close $TSV;
}

sub run_mummer_analysis {
    msg("running mummer analysis on blastn hits");
    
    my $jobnumber = 0;
    
    if (-s $suspect_reassign){
        unlink $suspect_reassign or err("failed to clean up temp file $suspect_reassign");
    }
    if (-s $mummer_log){
        unlink $mummer_log or err("failed to clean up $mummer_log");
    }
    
    CTG: foreach my $contig (sort(keys(%contigs))){
        next CTG if ($contigs{$contig}{REASSIGNED});
        next CTG if ($contigs{$contig}{ASSIGN});
        
        # skip if no hits
        if (!($contigs{$contig}{1})){
            $contigs{$contig}{ASSIGN} = "n";
            next CTG;
        }
        
        # skip if it's bigger than both it's hits
        if ($contigs{$contig}{2}){
            if ( $contigs{$contig}{LEN} > ( $contigs{$contigs{$contig}{1}}{LEN} + $contigs{$contigs{$contig}{2}}{LEN} ) ){
                $contigs{$contig}{ASSIGN} = 0;
                next CTG;
            }
        } elsif ( $contigs{$contig}{LEN} > $contigs{$contigs{$contig}{1}}{LEN} ){
            $contigs{$contig}{ASSIGN} = 0;
            next CTG;
        }
        
        # run the mummer job
        $available_threads->down(1);
        threads->create(\&mummer_job, $jobnumber, $contig);
        $jobnumber++;
    }
    
    # wait on remaining jobs
    sleep 1 while threads->list();
}
#-----
    sub mummer_job {
        my $tmp_log;
        my $cmd;
        
        my $job = $_[0];
        my $contig = $_[1];
        
        # make a ref seq
        my $ref = "$MUM_DIR/$job.ref.fa";
        qruncmd("cat $MINCE_DIR/$contigs{$contig}{1}.fasta > $ref");
        qruncmd("cat $MINCE_DIR/$contigs{$contig}{2}.fasta >> $ref") if ($contigs{$contig}{2});
        
        # run nucmer
        $cmd = "$Bin/../mummer/nucmer -p $MUM_DIR/$job.tmp $MINCE_DIR/$contig.fasta $ref 2>&1\n";
        $tmp_log .= "RUNNING: $cmd";
        $tmp_log .= `$cmd`;
        
        # run delta-filter
        $cmd = "$Bin/../mummer/delta-filter -m $MUM_DIR/$job.tmp.delta > $MUM_DIR/$job.tmp.m.delta\n";
        $tmp_log .= "RUNNING: $cmd";
        $tmp_log .= `$cmd`;
        
        $cmd = "$Bin/../mummer/delta-filter -r $MUM_DIR/$job.tmp.delta > $MUM_DIR/$job.tmp.r.delta\n";
        $tmp_log .= "RUNNING: $cmd";
        $tmp_log .= `$cmd`;
        
        # max-match coverage
        $cmd = "$Bin/../mummer/show-coords -b -c $MUM_DIR/$job.tmp.m.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'\n";
        $tmp_log .= $cmd;
        my $maxmatch = `$cmd`;
        $maxmatch =~ s/\s//g;
        $tmp_log .= "MAXMATCH coverage = $maxmatch\n";
        
        # best-align coverage
        $cmd = "$Bin/../mummer/show-coords -b -c $MUM_DIR/$job.tmp.r.delta | grep -P \"\\s+\\d\" | awk '{ s+=\$10 } END { print s }'\n";
        $tmp_log .= $cmd;
        my $bestmatch = `$cmd`;
        $bestmatch =~ s/\s//g;
        $tmp_log .= "BESTMATCH coverage = $bestmatch\n";
        
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
        
        # make a dotplot in unassigned
        if (-s "$UNASSIGNED/$contig.png"){
            unlink "$UNASSIGNED/$contig.png" or err("failed to clean up previous dotplot $UNASSIGNED/$contig.png");
        }
        if ($assignment ne "n"){
            $cmd = "$Bin/../mummer/mummerplot --fat -p $UNASSIGNED/$contig $MUM_DIR/$job.tmp.m.delta 2>&1\n";
            $tmp_log .= $cmd;
            $tmp_log .= `$cmd`;
            
            # cleanup
            foreach my $file("$UNASSIGNED/$contig.filter", "$UNASSIGNED/$contig.fplot", "$UNASSIGNED/$contig.rplot", "$UNASSIGNED/$contig.gp"){
                unlink $file or err("failed to clean up temp file $file\nERR_LOG:\n$tmp_log\n");
            }
        }
        
        #
        $writing_to_out->down(1);
        
        open my $TSV, ">>", $suspect_reassign or err("failed to open $suspect_reassign for appending");
        open my $MLOG, ">>", $mummer_log or err("failed to open $mummer_log for appending");
        
        if ($contigs{$contig}{2}){
            print $TSV "$contig\t$contigs{$contig}{1}\t$contigs{$contig}{2}\t$bestmatch\t$maxmatch\t$assignment\n";
        } else {
            print $TSV "$contig\t$contigs{$contig}{1}\t\t$bestmatch\t$maxmatch\t$assignment\n";
        }
        
        print $MLOG $tmp_log;
        
        close $TSV;
        close $MLOG;
        
        $writing_to_out->up(1);
        
        # clean up 
        foreach my $file ("$MUM_DIR/$job.ref.fa", "$MUM_DIR/$job.tmp.delta", "$MUM_DIR/$job.tmp.m.delta", "$MUM_DIR/$job.tmp.r.delta"){
            unlink $file or err("failed to clean up temp file $file");
        }
        
        # exit
        $available_threads->up(1);
        threads->detach();
    }
#-----


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
            
            # if on is repetitive, keep haplotig
            elsif ( $contigs{$ctg}{ASSIGN} . $contigs{$r_ctg}{ASSIGN} =~ /rh|hr/i){
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
        if ($junk{$ctg}){
            qruncmd("cat $MINCE_DIR/$ctg.fasta >> $out_artefacts");
        } elsif ($contigs{$ctg}{REASSIGNED}){
            write_seq($CUH, $ctg);
            if ($contigs{$ctg}{2}){
                print $CUT "$ctg\t$contigs{$ctg}{1}\t$contigs{$ctg}{2}\t$contigs{$ctg}{BM}\t$contigs{$ctg}{MM}\t$contigs{$ctg}{ASSIGN}\n";
            } else {
                print $CUT "$ctg\t$contigs{$ctg}{1}\t.\t$contigs{$ctg}{BM}\t$contigs{$ctg}{MM}\t$contigs{$ctg}{ASSIGN}\n";
            }
        } else {
            qruncmd("cat $MINCE_DIR/$ctg.fasta >> $out_fasta");
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


#---UTILITY SUBROUTINES---

sub print_message {
    my $t = localtime;
    my $line = $t->dmy . " " . $t->hms . " @_\n";
    print STDERR $line;
    print $LOG $line if ($LOG);
}

sub msg {
    print_message("INFO: @_");
}

sub err {
    print_message("ERROR: @_\n\nPurge_Haplotigs has failed.\n");
    exit(1);
}

sub runcmd {
    print_message("RUNNING: @_");
    system("@_") == 0 or err("Failed to run @_\nCheck $_[-1] for possible causes.");
    print_message("FINISHED: @_")
}

sub qruncmd {
    system("@_") == 0 or err("Failed to run @_");
}

sub check_files {
    my $check=1;
    foreach(@_){
        if (!(-s $_)){
            print_message("ERROR: file \"$_\" does not exist or is empty");
            $check=0;
        }
    }
    return $check;
}

sub chkprog {
    print_message("CHECKING $_[0]");
    system("@_") == 0 or (print_message("ERROR: missing program $_[0]"), return 0);
    msg("$_[0] OK");
    return 1;
}



