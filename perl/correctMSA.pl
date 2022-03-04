use FileUtil;
use BioStats;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

correctMSA.pl - Take an MSA and identify and correct potential errors

=head1 VERSION

This document refers to version 1.00 of correctMSA.pl, released 12.27.2021.

=head1 SYNOPSIS

correctMSA.pl [-help]

Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message
    -reset                Overwrite existing results
    -noN                  Do not replace bases with N's during correction

    -msa                  Input MSA file
    -out                  Write output file here
    -rocitsPath           Full path to where RoC-ITS scripts are located

    -minCov               Minimum coverage threshold that distinguishes noise from a meaningful clade (default = 0)
    -maxEdits             Maximum number of edits allowed before read is rejected (default 0.5% of average length)

    -reference            (Optional) Reference indicating safe identifiers

=head1 DESCRIPTION

=head2 Overview

correctMSA.pl is designed to correct problems with an improved MSA that arise due to errors
in the RoC-ITS read consensus. Minor indels and substitutions occur in the 1% range and this
tool attempts to identify they "rare" errors and when possible "fix" them towards the consensus
to reduce noise in the data.

In some cases, the MSA may include reference sequences. In these cases it would be
inappropriate to "fix" a known reference sequence. These sequences can be specified by a 
regular expression and they will not be used during the "correction" process.

The minCov parameter is a threshold that distinguishes meaningful data from suspected noise/error.
Setting this value is crucial. By default, all the data is kept. Any value below 2 will keep all
the data and will ensure that all the noise/error is shown. Values equal to or greater than 2 will
start to remove noise/error. Here is a practical example. Given 100 reads from an isolate with 10
ribosomal operons, its safe to assume that each operon should have 5 or more reads. So setting a
minCov of 5 would eliminate much of the noise without losing any individual rrns. However, this may
be too stringent for a metagenomic sample where there could be two or more distinct organisms
present so this approach should be used with care to avoid "correcting" away real differences
in the data.

The maxEdits removes sequences from the output that have more edits than expected. This tends to
remove sequences with large indels/chimera. The log will report how many sequences were excluded
using this parameter which could also be used to identify clusters consisting of multiple distinct
organisms. The default parameter was set based on usage and empirical data from known isolates.

=cut

my ($help, $log, $msa, $output, $minCov, $excluded, $maxEdits, $reference, $noNs);
my ($rocitsPath);

BEGIN: {
  GetOptions(
    "help|h"        => \$help,
    "log"           => \$log,
    "msa=s"         => \$msa,
    "out=s"         => \$output,
    "excluded=s"    => \$excluded,
    "minCov=s"      => \$minCov,
    "noN"           => \$noNs,
    "rocitsPath=s"  => \$rocitsPath,
    "maxEdits:s"    => \$maxEdits,
    "reference:s"   => \$reference,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}
if (!defined($msa) || !defined($output) || !defined($excluded)) {
  pod2usage(1);
  exit;
}

my $gggs = "n" x 22;
$reference = "Unknown" unless defined($reference);

### read in MSA data
my @lens;
my $mfh = FileUtil::openFileHandle("perl $rocitsPath/reformatFasta.pl $msa 0 99999999|");
while (<$mfh>) {
  if (/^>(\S+)/) {
    my $id = $1;
    $_ = <$mfh>;
    chomp;
    my $s = $_;
    $s =~ s/-//g;
    $s =~ s/n//g;
    push @lens,length($s);
    $seq{$id} = $gggs.uc($_).$gggs;
    my @seq = split //,$seq{$id};
    $len = @seq;
    push @seqs1,\@seq;
    push @ids,$id;
  }
}
$mfh->close();

my $avgLen = BioStats::getMean(\@lens);
$maxEdits = $avgLen * 0.005 unless defined($maxEdits);
print STDERR "LOG: Max edits allowed = ",$maxEdits,": $0\n" if $log;

### look at the number of bases in a particular column of the MSA and reject that column if there are differences but
### the number of differences is less than minCov (includes gaps)
my @good;
my $alteredBases = 0;
my $cnt = @ids;
for (my $i = 0; $i < $len; $i++) {
  my $bases = 0;
  my %dist;
  for (my $j = 0; $j < $cnt; $j++) {
    $dist{${$seqs1[$j]}[$i]}++;
  }
  my $best;
  ### how many different bases were seen here
  my $num = keys %dist;
  if ($num == 1) {
    ### if there was only a single base in this column...
    $bases = 999999999;
  } else {
    ### if there was more than one base sort the bases from most to least abundant
    my $first = 0;
    foreach my $k (sort {$dist{$b} <=>$dist{$a}} keys %dist) {
      if ($first) {
        ### count all the non-most abundant bases
        $bases+= $dist{$k};
      } else {
        ### the most abundant base is tagged
        $first++;
        $best = $k;
      }
    }
  }
  ### Here editting takes place. Editting should be conservative and edits are meant to deal with
  ### noise/error in the read data. The minCov crucially sets this threshold.
#  print "xx $i $best $bases $minCov\n";
  
  
  ### check if the number of non-best bases is sufficiently large to be worth noting
  if ($bases < $minCov) {
    ### if its below the minCov then just replace all the non-canonical bases with the best
    for (my $j = 0; $j < $cnt; $j++) {
      if (${$seqs1[$j]}[$i] ne $best) {
        if ($ids[$j] =~ /$reference/o) {
          ### don't alter the references IDs just lowercase differences
          ${$seqs1[$j]}[$i] = lc(${$seqs1[$j]}[$i]);
        } else {
          ${$seqs1[$j]}[$i] = lc($best);
          $altered{$ids[$j]}++;
          $alteredBases++;
        }
      }
    }
  } else {
    ### if the number of not best bases is sufficient decide whether each base is above the threshold
    ### if the number is above the minCov threshold keep the base whatever it is
    ### if not above the minCov threshold then we know there are 2 or more possibilities so replace
    ###   base with an N to indicate that we think there is an error here but we can't correct it
    for (my $j = 0; $j < $cnt; $j++) {
#  print "yy $i $j $best ${$seqs1[$j]}[$i] $bases $minCov $dist{${$seqs1[$j]}[$i]}\n";
      if ($dist{${$seqs1[$j]}[$i]} >= $minCov) {
      } else {
        if ($ids[$j] =~ /$reference/o) {
        } else {
          $altered{$ids[$j]}++;
          if ($noNs) {
          } else {
            ${$seqs1[$j]}[$i] = "N";
          }
          $alteredBases++;
        }
      }
    }
  }
}
print STDERR "LOG: Altered $alteredBases non-reference bases: $0 line ",__LINE__,"\n" if $log;

### write output
my $editedOut = 0;
my $fho = new FileHandle(">$output");
my $fho2 = new FileHandle(">$excluded");
for (my $j = 0; $j < $cnt; $j++) {
  my $tmp = join "",@{$seqs1[$j]};
  my $n = $altered{$ids[$j]}+0;
  if ($n <= $maxEdits) {
    print $fho ">$ids[$j] fixed=$n\n$tmp\n";
  } else {
    print $fho2 ">$ids[$j] fixed=$n $maxEdits x\n$tmp\n";
    $editedOut++;
  }
  
}
$fho->close();
$fho2->close();
print STDERR "LOG: MaxEdit removed $editedOut sequences of $cnt sequences: $0 line ",__LINE__,"\n" if $log;
exit;
