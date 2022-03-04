use FileUtil;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

improveMSA.pl - Takes a set of genus specific reads, cluster them at a specific percent
  identity cutoff, and generate an MSA and consensus sequence for each cluster.

=head1 VERSION

This document refers to version 1.00 of improveMSA.pl, released 12.27.2021.

=head1 SYNOPSIS

improveMSA.pl [-help]

Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message
    -reset                Overwrite existing results

    -msa                  Input MSA file
    -out                  Write output file here
    -rocitsPath           Full path to where RoC-ITS scripts are located

    -minimumCov           Minimum coverage required to keep a column in the MSA (default 0)
    -wordSize             Anchor size for conserved regions (default 12)
    -maxJobs              Max number of concurrent jobs (default 40)
    -local                Directory where temporary files will be written (default ".")

    -probcons             (Optional) Full path to probcons executable

=head1 DESCRIPTION

=head2 Overview

improveMSA.pl is designed to correct problems with muscle MSA. The probcons MSA tool tends to
correct these problems but is MUCH slower. Here we identify anchor regions (regions that are
nearly identical in the vast majority of reads) and then re-align the intervening bases. Anchors
must have a reasonable minimum lenght to be effective (8-20 bps). The larger the word size, the
fewer blocks of sequence there will be, and they will be appropriately larger slowing down the
re-alignment process. Too small and the blocks are not guaranteed to be meaningful and there is
a possibility that the alignment will be sub-optimal.

=cut

my ($help, $log, $input, $output, $word, $minCov, $jobs, $localDir);
my ($probconsBin, $rocitsPath);

BEGIN: {
  GetOptions(
    "help|h"        => \$help,
    "log"           => \$log,
    "msa=s"         => \$input,
    "out=s"         => \$output,
    "wordSize:s"    => \$word,
    "local:s"       => \$localDir,
    "rocitsPath=s"  => \$rocitsPath,
    "maxJobs:s"     => \$jobs,
    "minimumCov:s"  => \$minCov,
    "probcons:s"    => \$probconsBin,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}
if (!defined($input)) {
  pod2usage(1);
  exit;
}

$word = 12 unless defined($word);
$minCov = 0 unless defined($minCov);
$localDir = "." unless defined($localDir);
$probconsBin = "probcons" unless defined($probconsBin);
$jobs = 40 unless defined($jobs);
my $gggs = "n" x 22;

my $mfh = FileUtil::openFileHandle("perl $rocitsPath/reformatFasta.pl $input 0 99999999|");
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
  my $bases = 0;
  for (my $j = 0; $j < $cnt; $j++) {
    $dist{${$seqs1[$j]}[$i]}++;
    if (${$seqs1[$j]}[$i] eq "-") {
    } else {
      $bases++;
    }
  }
  if ($bases < $minCov) {
    for (my $j = 0; $j < $cnt; $j++) {
      ${$seqs1[$j]}[$i] = "";
    }
  }
}
for (my $j = 0; $j < $cnt; $j++) {
  my $tmp = join "",@{$seqs1[$j]};
  my @seq = split //,$tmp;
  $len = @seq;
  push @seqs,\@seq;
}

for (my $i = 0; $i < $len; $i++) {
  my %dist;
  for (my $j = 0; $j < $cnt; $j++) {
    $dist{${$seqs[$j]}[$i]}++;
  }
  if ($dist{"-"}/$cnt >= 0.1) {
    push @good,"-";
  } else {
    my $best;
    foreach my $j (sort {$dist{$b}<=>$dist{$a}; } keys %dist) {
      $best = $j;
      last;
    }
    if ($dist{$best}/$cnt < 0.9) {
      push @good,"-";
    } else {
      push @good,"G";
    }
  }
}

my $g = join "",@good;
if ($g =~ /^(-+)/) {
  my $w = "p" x length($1);
  $g =~ s/^(-+)/$w/;
}

my @gg = split //,$g;
my @regionBegins;
my @regionEnds;
my $start = 0;
my $endFound = 0;
while (!$endFound) {
  my $i = $start;
  for ($i = $start; $i < @gg; $i++) {
    if ($gg[$i] eq "-") {
      $begin = $i - $word;
      my $tot = 0;
      my $pos = $i;
      while ($tot < $word) {
        if ($gg[$pos] eq "G") {
          $tot++;
        } else {
          $tot = 0;
        }
        $pos++;
      }
      $end = $pos;
#      my $x = join "",@gg[$begin..$end];
#      print "$i $begin $end\n$x\n";
      push @regionBegins,$begin;
      push @regionEnds,$end;
      $start = $end - $word;
      last;
    }
  }
  if ($i < @gg) {
  } else {
    $endFound = 1;
  }
}

for (my $ii = 0; $ii < @regionBegins; $ii++) {
  my $fho = new FileHandle(">$localDir/temp.$$.$ii.$regionBegins[$ii].$regionEnds[$ii].fa");
  for (my $i = 0; $i < $cnt; $i++) {
    my @set;
    for (my $j = $regionBegins[$ii]; $j < $regionEnds[$ii]; $j++) {
      my $chr = ${$seqs[$i]}[$j];
      if ($j < ($regionEnds[$ii]-$word)) {
        ${$seqs[$i]}[$j] = "";
        push @set,$chr;
      } else {
        ### subsequences that have gaps in the right anchor region are by definition rare but need to
        ### be replaced by N's. This is important as later non-dashes in the right anchor are removed
        ### to when the new msa is inserted.
        if ($chr eq "-") {
          push @set,"n";
        } else {
          push @set,$chr;
        }
      }
    }
    my $sub = join "",@set;
    $sub =~ s/-+//g;
    print $fho ">$ids[$i] $regionBegins[$ii] $regionEnds[$ii]\n";
    print $fho $sub,"\n";
  }
  $fho->close();
  my $probCmd = "$probconsBin -c 0 $localDir/temp.$$.$ii.$regionBegins[$ii].$regionEnds[$ii].fa |perl $rocitsPath/reformatFasta.pl - 0 9999999 > $localDir/temp.$$.$ii.$regionBegins[$ii].out &";
  push @cmds,$probCmd;
  print STDERR $probCmd,"\n" if $log;
}
my $ccnt = 0;
my $fho = new FileHandle(">$localDir/temp.$$.cmd");
foreach my $i (@cmds) {
  print $fho $i,"\n";
  $ccnt++;
  if ($cnt == $jobs) {
    print $fho "wait\n";
    $ccnt = 0;
  }
}
print $fho "wait\n";
$fho->close();

{
  ### have to go to some lengths to redirect stderr not to clutter output
  my $cmd = "tcsh -x $localDir/temp.$$.cmd > /dev/null";
  print STDERR $cmd,"\n";
  open(my $STDOLD, '>&', STDERR);
  open(STDERR, '>>', "/dev/null");
  my $run = system($cmd);
  open(STDERR, '>&', $STDOLD);
}

for (my $ii = 0; $ii < @regionBegins; $ii++) {
  my @result;
  my $rfh = FileUtil::openFileHandle("$localDir/temp.$$.$ii.$regionBegins[$ii].out");
  while (<$rfh>) {
    next if /^>/;
    chomp;
    push @result,$_;
  }
  $rfh->close();
  for (my $i = 0; $i < $cnt; $i++) {
    for (my $j = $regionBegins[$ii]; $j < $regionEnds[$ii]; $j++) {
      my @res = split //,$result[$i];
      my $e = @res;
      my $removed = 0;
      for (my $k = $e - 1; $k >= 0; $k--) {
        if ($res[$k] =~ /\w/o) {
          $res[$k] = "";
          $removed++;
        }
        last if $removed == $word;
      }
      my $res = join "",@res;
      ${$seqs[$i]}[$j] = lc($res);
      last;
    }
  }
}

my $fho = new FileHandle(">$output");
for (my $i = 0; $i < $cnt; $i++) {
  my $l = join "",@{$seqs[$i]};
  $l =~ s/n/-/g;
  print $fho ">$ids[$i]\n$l\n";
}
$fho->close();

$cmd = "rm $localDir/temp.$$*.fa";
system($cmd);
$cmd = "rm $localDir/temp.$$*.out $localDir/temp.$$.cmd";
system($cmd);

exit;
