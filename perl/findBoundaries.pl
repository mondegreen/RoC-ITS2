use FileUtil;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

findBoundaries.pl - takes an MSA and identifies poorly represented regions on the ends of the MSA that
can be trimmed.

=head1 VERSION

This document refers to version 1.00 of processRoCITSSampleSet.pl, released 12.27.2021.

=head1 SYNOPSIS

findBoundaries.pl [-help]

Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message

    -msa                  Input fasta format MSA file
    -covCutoff            Number between 0 and 1 inclusive; default 0.98
    -rocitsPath           Full path to where RoC-ITS scripts are located

=head1 DESCRIPTION

=head2 Overview

The covCutoff parameter specifies the stopping condition for determining the boundaries. The value,
between 0 and 1 inclusive, specifies the percentage of usable bases required to determine the boundary.
A usable base is an A, T, C or G. Unusable bases include gaps or N's.

It prints the length of the two bounds as a space separated tuple.

=cut

my ($help, $log, $fasta, $per);

BEGIN: {
  GetOptions(
    "help|h"        => \$help,
    "log"           => \$log,
    "msa=s"         => \$fasta,
    "covCutoff:s"   => \$per,
    "rocitsPath=s"  => \$rocitsPath,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}

$per = 0.98 unless defined($per);

my $mfh = FileUtil::openFileHandle("perl /N/home/drusch/jcviBin/altBin/reformatFasta.pl $fasta 0 99999999|");
while (<$mfh>) {
  if (/^>(\S+)/) {
    my $id = $1;
    $_ = <$mfh>;
    chomp; 
    $seq{$id} = $gggs.uc($_).$gggs;
    my @seq = split //,$seq{$id};
    $len = @seq;
    push @seqs1,\@seq;
    push @ids,$id;
  }
}
$mfh->close();

my @good;
my $cnt = @ids;
for (my $i = 0; $i < $len; $i++) {
  my $good = 0;
  my $bad = 0;
  for (my $j = 0; $j < $cnt; $j++) {
    my $b = ${$seqs1[$j]}[$i];
    if ($b eq "N" or $b eq "-" or $b eq "n") {
      $bad++;
    } else {
      $good++;
    }
  }
#  print "$i\t$good\t$cnt\t$per\n";
  if ($good/$cnt < $per) {
    $good[$i] = 0;
  } else {
    $good[$i] = 1;
  }
}
my $good = join "",@good;
#print $good,"\n";
if ($good =~ /^(0*)/) {
  $b = $1;
}
if ($good =~ /.*1(0*)$/) {
  $e = $1;
}
my ($bb,$ee,$bbl,$eel);
($bbl,$eel) = (9999,9999);
if (defined($b) && defined($e)) {
  my $bl = length($b);
  my $el = length($e);
  if ($bl < $bbl) {
    $bbl = $bl;
    $bb = $b;
  }
  if ($el < $eel) {
    $eel = $el;
    $ee = $e;
  }
}
$bb =~ tr/0/N/;
$ee =~ tr/0/N/;

print join(" ",length($bb),length($ee)),"\n";
