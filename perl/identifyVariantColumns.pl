use FileUtil;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

identifyVariantColumns.pl - Takes an MSA and looks at the variation in each column. When the variation is
equal to or above a percentage cutoff, then it marks that as a column of interest. The columns of interest
are then extracted to make subsequences that can then be clustered.

=head1 VERSION

This document refers to version 1.00 of correctMSA.pl, released 12.27.2021.

=head1 SYNOPSIS

identifyVariantColumns.pl [-help]

Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message

    -msa                  Input MSA file
    -out                  Write output file here
    -leftBound            Trim N bases from 5' end of MSA
    -rightBound           Trim N bases from 3' end of MSA
    -rocitsPath           Full path to where RoC-ITS scripts are located

    -pid                  Cutoff to decide if a column is informative
    -context              Number of bases to select around column of interest (default = 2)

=head1 DESCRIPTION

=head2 Overview

identifyVariantColumns.pl

The number of variable positions in RoC-ITS reads within a given species can be very low (fractions of
a percent). By identifying only the variant columns of the MSA and using these to construct subsequences
for clustering, less nuanced cutoffs can be used to identify meaningful clusters in the data.

The MSA should be trimmed otherwise the ends of the MSA disproportionately dominate the clustering.

Contextual information can help to balance and expand the number of bases that are extracted.

If the number of identified columns is below a specified cutoff, then the entire MSA is marked.

The subsequences produced are no longer aligned as gaps are not included.

=cut

my ($help, $log, $msa, $output, $ignoreBegin, $ignoreEnd, $perc, $minColumns, $context);
my ($rocitsPath);

BEGIN: {
  GetOptions(
    "help|h"        => \$help,
    "log"           => \$log,
    "msa=s"         => \$msa,
    "output=s"      => \$output,
    "pid=s"         => \$perc,
    "leftTrim=s"    => \$ignoreBegin,
    "rightTrim=s"   => \$ignoreEnd,
    "minColumns=s"  => \$minColumns,
    "rocitsPath=s"  => \$rocitsPath,
    "context:s"     => \$context,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}

$context = 2 unless defined($context);

### dump indicates whether to output the subsequences, value indicates how many bases to offset for the
### perc value provided. Not that negative values shift towards the 5' end while positive values shift 
### towards the 3' end

my $ffh = FileUtil::openFileHandle("perl $rocitsPath/reformatFasta.pl $msa 0 999999999 1|");
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
my $intCount = 0;
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
    my $per = $b{$i}/$num;
    if ($per < $perc) {
      $interesting[$j] = 1;
      $intCount++;
    } else {
      $interesting[$j] = 0;
    }
    last;
	}
}
print STDERR "LOG: Interesting columns = ",$intCount," out of $len: $0 line ",__LINE__,"\n" if $log;


my @interesting2;
for (my $i = 0; $i < $len; $i++) {
  if ($interesting[$i]) {
    my $max = $i + $context;
    for (my $j = $i - $context; $j <= $max; $j++) {
      if ($j >= $ignoreBegin && $j < $len-$ignoreEnd) {
        $interesting2[$j] = 1;
      }
    }
  }
}

if ($intCount < $minColumns) {
  print STDERR "LOG: Not enough interesting columns detected. Using entire sequence: $0 line ",__LINE__,"\n";
  for (my $j = $ignoreBegin; $j < $len-$ignoreEnd; $j++) {
    $interesting2[$j] = 1;
  }
}

my $fho = new FileHandle(">$output");
for (my $i = 0; $i < @ids; $i++) {
  print $fho ">$ids[$i]\n";
  my @s;
  for (my $j = $ignoreBegin; $j < ($len-$ignoreEnd); $j++) {
    push @s,${$seq[$i]}[$j] if $interesting2[$j];
  }
  my $l = join "",@s;
  $l =~ s/-//g;
  print $fho $l,"\n";
}
$fho->close();
