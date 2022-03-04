use FileUtil;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

extractEffectiveColumns.pl - Takes an MSA, trims the uneven edges, and produces a fasta file of trimmed fasta

=head1 VERSION

This document refers to version 1.00 of correctMSA.pl, released 12.27.2021.

=head1 SYNOPSIS

extractEffectiveColumns.pl [-help]

Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message

    -msa                  Input MSA file
    -out                  Write output file here
    -leftBound            Trim N bases from 5' end of MSA
    -rightBound           Trim N bases from 3' end of MSA
    -rocitsPath           Full path to where RoC-ITS scripts are located

=head1 DESCRIPTION

=head2 Overview

extractEffectiveColumns.pl

The effective columns are all the columns of the MSA after trimming. This script performs the trimming
and prints out the resulting fasta subsequences.

If the number of identified columns is below a specified cutoff, then the entire MSA is marked.

The subsequences produced are no longer aligned as gaps are not included.

=cut

my ($help, $log, $msa, $output, $ignoreBegin, $ignoreEnd, $perc, $minColumns);
my ($rocitsPath);

BEGIN: {
  GetOptions(
    "help|h"        => \$help,
    "log"           => \$log,
    "msa=s"         => \$fasta,
    "output=s"      => \$output,
    "leftTrim=s"    => \$ignoreBegin,
    "rightTrim=s"   => \$ignoreEnd,
    "rocitsPath=s"  => \$rocitsPath,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}

$perc = 2;

my $ffh = FileUtil::openFileHandle("perl $rocitsPath/reformatFasta.pl $fasta 0 999999999 1|");
while (<$ffh>) {
	if (/^>(\S+)/) {
		$id = $1;
		$cnt++;
	} else {
		chomp;
		my $ll = length($_);
		my ($dash) = $_ =~ tr/-/-/;
		$totalChar += $ll;
		my @s = split //;
		my $l = @s;
		$len = $l if $l > $len;
		push @seq,\@s;
		push @ids,$id;
	}
}
$ffh->close();

my $totalMismatches = 0;
my $totalPositions = 0;
my $num = @ids;
my @con;
for (my $j = $ignoreBegin; $j < $len-$ignoreEnd; $j++) {
	my %b;
  my $n = 0;
	for (my $i = 0; $i < @ids; $i++) {
		$b{${$seq[$i]}[$j]}++;
    $n++;
    $totalPositions++;
	}
  $con[$j] = "-";
	foreach my $i (sort { $b{$b} <=> $b{$a};} keys %b) {
    if (($b{$i}/$num) < $perc) {
      $interesting[$j] = 1;
    }
    last;
	}
}

my $fho = new FileHandle(">$output");
for (my $i = 0; $i < @ids; $i++) {
  print $fho ">$ids[$i]\n";
  my @s;
  for (my $j = $ignoreBegin; $j < ($len-$ignoreEnd); $j++) {
    push @s,${$seq[$i]}[$j] if $interesting[$j];
  }
  my $l = join "",@s;
  $l =~ s/-//g;
  print $fho $l,"\n";
}
$fho->close();
