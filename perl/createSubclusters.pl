use FileUtil;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

createSubclusters.pl - given a fasta file and clusters given by a modified cd-hit output,
create a set of fasta files for each cluster found.

=head1 VERSION

This document refers to version 1.00 of correctMSA.pl, released 12.27.2021.

=head1 SYNOPSIS

correctMSA.pl [-help]

Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message
    -reset                Overwrite existing results

    -selected             Clusters that have been selected for keeping
    -clstr                Input cd-hit .clstr file
    -fasta                Fasta file to be distributed into clusters
    -baseName             Base filename that the clusters will be written to
    -addCount             Add a tag to the cluster file name with the number of sequences included (default = 0)
    -centroidOnly         Only output the reference sequence identified by cd-hit with a '*'
    -rocitsPath           Full path to where RoC-ITS scripts are located

=head1 DESCRIPTION

=head2 Overview

=cut

my ($help, $log, $cls, $selected, $clstr, $fasta, $base, $noCnt, $centroidOnly, $rocitsPath);
my ($rocitsPath);

BEGIN: {
  GetOptions(
    "help|h"        => \$help,
    "log"           => \$log,
    "selected=s"    => \$pick,
    "clstr=s"       => \$cls,
    "fasta=s"       => \$fasta,
    "baseName=s"    => \$base,
    "addCount:s"    => \$noCnt,
    "centroidOnly:s"=> \$centroidOnly,
    "rocitsPath=s"  => \$rocitsPath,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}

my %clstrs;
my $pfh = FileUtil::openFileHandle($pick);
while (<$pfh>) {
  chomp;
  if (/^(\d+)/) {
    $clstrs{$1} = 1;
  }
}
$pfh->close();

my $cfh = FileUtil::openFileHandle($cls);
while (<$cfh>) {
  if (/^>Cluster\s+(\d+)/) {
    $c = $1;
  } else {
    /^\S+\s+\S+\s+>(\S+)\.\.\./;
    my $id = $1;
    if ($centroidOnly) {
      if (/\.\.\. \*/) {
        $m{$c}{$id} = 1;
      }
    } else {
      $m{$c}{$id} = 1;
    }
  }
}
$cfh->close();

my $ffh = FileUtil::openFileHandle($fasta);
while (<$ffh>) {
  if (/^>(\S+)/) {
    my $id = $1;
    $_ = <$ffh>;
    chomp;
    $seq{$id} = $_;
  }
}
$ffh->close();

foreach my $i (keys %clstrs) {
  my $cnt = keys %{$m{$i}};
  my $fout;
  if (defined($noCnt) && $noCnt) {
    $fout = new FileHandle(">$base.$i.fa");
  } else {
    $fout = new FileHandle(">$base.$i.num=$cnt.fa");
  }
  foreach my $j (keys %{$m{$i}}) {
    print $fout ">$j\n$seq{$j}\n";
  }
  $fout->close();
}
