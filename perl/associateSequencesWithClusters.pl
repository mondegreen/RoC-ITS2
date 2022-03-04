use FileUtil;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

associateSequencesWithClusters.pl - reformat the output of cd-hit clustering and associate each
sequence ID with a particular cluster while filtering out small clusters

=head1 VERSION

This document refers to version 1.00 of correctMSA.pl, released 12.27.2021.

=head1 SYNOPSIS

correctMSA.pl [-help]

Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message
    -reset                Overwrite existing results

    -clstr                Input cd-hit .clstr file
    -rocitsPath           Full path to where RoC-ITS scripts are located
    -minSize              Minimum size of a cluster (default = 0)
    
    -reference            Sequences to be saved regardless of cluster size
    -keyword

=head1 DESCRIPTION

=head2 Overview

=cut

my ($help, $log, $f, $mmin, $keyword, $rocitsPath, $reference);
my ($rocitsPath);

BEGIN: {
  GetOptions(
    "help|h"        => \$help,
    "log"           => \$log,
    "clstr=s"       => \$f,
    "minSize=s"     => \$mmin,
    "rocitsPath=s"  => \$rocitsPath,
    "reference:s"   => \$reference,
    "keyword:s"     => \$keyword,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}

my $fh = FileUtil::openFileHandle($f);
while (<$fh>) {
  if (/^>Cluster\s+(\d+)/) {
    $c = $1;
  } else {
    $cnt{$c}++;
    $total++;
  }
}
$fh->close();

my $min = 99999999;
my $fh = FileUtil::openFileHandle($f);
while (<$fh>) {
  if (/^>Cluster\s+(\d+)/) {
    $c = $1;
  } else {
    if (/$keyword/o) {
      $min = $cnt{$c} if $cnt{$c} < $min;
    }
  }
}
$fh->close();

if ($mmin && !defined($keyword)) {
  $min = $mmin;
}
print STDERR "LOG: Minimum cluster used = ",$min,": $0\n" if $log;
my $fh = FileUtil::openFileHandle($f);
while (<$fh>) {
  if (/^>Cluster\s+(\d+)/) {
    $c = $1;
    if ($cnt{$c} >= $min) {
      $seen += $cnt{$c};
    }
  } else {
    if ($cnt{$c} >= $min) {
      print "$c\t$cnt{$c}\t$_";
    } else {
      if (/$reference/o) {
        print "$c\t$cnt{$c}\t$_";
      }
    }
  }
}
$fh->close();

my $per = int(1000*$seen/$total)/10;
print STDERR "LOG: $seen out of $total sequences kept ($per%) = ",$min,": $0\n" if $log;
