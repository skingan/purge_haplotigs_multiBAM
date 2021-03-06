package PipeUtils;

use strict;
use warnings;
use Time::Piece;

use Exporter 'import';

our @EXPORT = qw(msg err runcmd qruncmd check_programs check_files);


sub print_message {
    my $t = localtime;
    my $line = $t->dmy . " " . $t->hms . " @_\n";
    print STDERR $line;
    print $::LOG $line if ($::LOG);
}

sub msg {
    print_message("INFO: @_");
}

sub err {
    print_message("ERROR: @_\n\nPIPELINE FAILURE\n");
    exit(1);
}

sub runcmd {
    my $job = shift;
    ($job->{silent}) || print_message("RUNNING: $job->{command}");
    if (system("$job->{command}") != 0){
        print_message("ERROR: Failed to run $job->{command}");
        print_message("Check $job->{logfile} for possible errors") if ($job->{logfile});
        err("Exiting due to job failure");
    } else {
        ($job->{silent}) || print_message("FINISHED: $job->{command}");
    }
}

sub qruncmd {
    system(@_) == 0 or err("Failed to run @_\n");
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

sub check_programs {
    my $chk=1;
    my $t = localtime;
    my $line = $t->dmy . " " . $t->hms;
    foreach my $prog (@_){
        print STDERR "$line CHECKING $prog... ";
        my $notexists = `type $prog 2>&1 1>/dev/null || echo 1`;
        if ($notexists){
            print STDERR "ERROR: missing program $prog\n";
            $chk = 0;
        } else {
            print STDERR "OK\n";
        }
    }
    return $chk;
}
