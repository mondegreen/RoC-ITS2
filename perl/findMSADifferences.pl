use FileUtil;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

findMSADifferences.pl - Identify salient differences in an MSA

=head1 VERSION

This document refers to version 1.00 of processRoCITSSampleSet.pl, released 12.27.2021.

=head1 SYNOPSIS

findMSADifferences.pl [-help]

Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message

    -msa                  Input fasta format MSA file
    -covCutoff            Clustering percent identity (times 10 so 970 equals 97.0%)
    -rocitsPath           Full path to where RoC-ITS scripts are located
    
    -interesting          Only show columns that are above the coverage cutoff
    -ignoreCols           Optional comma separated list; Ignore these columns by replacing with - 

=head1 DESCRIPTION

=head2 Overview

The covCutoff parameter specifies the stopping condition for determining the boundaries. The value,
between 0 and 1 inclusive, specifies the percentage of usable bases required to determine the boundary.
A usable base is an A, T, C or G. Unusable bases include gaps or N's.

It prints the length of the two bounds as a space separated tuple.

=cut

my ($help, $log, $msa, $per, $showInterestingOnly, $ignored, $showInt, $showPos);

BEGIN: {
  GetOptions(
    "help|h"        => \$help,
    "log"           => \$log,
    "msa=s"         => \$msa,
    "covCutoff:s"   => \$min,
    "interesting"   => \$showInterestingOnly,
    "ignoreCols:s"  => \$ignored,
    "showInteresting" => \$showInt,
    "showPos"       => \$showPos,
    "rocitsPath=s"  => \$rocitsPath,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}

$min = 0 unless defined($min);
$showInterestingOnly = 0 unless defined($showInterestingOnly);
my @ignored = split /,/,$ignored;

my $mfh = FileUtil::openFileHandle("perl $rocitsPath/reformatFasta.pl $msa 0 99999999|");
while (<$mfh>) {
  if (/^>(\S+)/) {
    my $id = $1;
    $def{$id} = $_;
    $_ = <$mfh>;
    chomp; 
    $seq{$id} = uc($_);
    my @seq = split //,$seq{$id};
    foreach my $i (@ignored) {
      $seq[$i] = "-";
    }
    $seq{$id} = join "",@seq;
    if (!$showInterestingOnly) {
      print $def{$id};
      print $seq{$id},"\n";
    }
    $len = @seq;
    push @seqs,\@seq;
    push @ids,$id;
  }
}
$mfh->close();

my @good;
my $cnt = @ids;
for (my $i = 0; $i < $len; $i++) {
  my %dist;
  for (my $j = 0; $j < $cnt; $j++) {
    $dist{${$seqs[$j]}[$i]}++;
  }
  my $best;
  foreach my $j (sort {$dist{$b}<=>$dist{$a}; } keys %dist) {
    $best = $j;
    last;
  }
  if ($dist{$best}/$cnt < (1-$min)) {
    push @good,"-";
  } else {
    push @good,"G";
  }
}

my %pos;
if ($showInterestingOnly) {
  for (my $j = 0; $j < $cnt; $j++) {
    my @s;
    for (my $i = 0; $i < $len; $i++) {
      push @s,${$seqs[$j]}[$i] if ($good[$i] eq "-");
    }
    my $l = join "",@s;
    print ">$ids[$j]\n$l\n";
  }
  my @s;
  for (my $i = 0; $i < $len; $i++) {
    if ($good[$i] eq "-") {
      push @s,$good[$i];
      $pos{$i} = 1;
    }
  }
  if ($showInt) {
    my $l = join "",@s;
    print "$l\n";
  }
  my $l = join ",",@p;
  print $l,"\n";
} else {
  if ($showInt) {
    my $g = join "",@good;
    print $g,"\n";
  }
}
if ($showPos) {
  my @row;
  for (my $i = 0; $i < $len; $i++) {
    my $n = sprintf "%4s",$i;
    my @num = split //,$n;
    if ($showInterestingOnly) {
      if ($good[$i] eq "-") {
        push @{$row[0]},$num[0]; 
        push @{$row[1]},$num[1]; 
        push @{$row[2]},$num[2]; 
        push @{$row[3]},$num[3]; 
      }
    } else {
      push @{$row[0]},$num[0]; 
      push @{$row[1]},$num[1]; 
      push @{$row[2]},$num[2]; 
      push @{$row[3]},$num[3]; 
    }
  }
  for (my $i = 0; $i < 4; $i++) {
    my $l = join "",@{$row[$i]};
    my $x = $l =~ tr/ / /;
    if ($x == length($l)) {
      ### don't print all spaces
    } else {
      print $l,"\n";
    }
  }
}
