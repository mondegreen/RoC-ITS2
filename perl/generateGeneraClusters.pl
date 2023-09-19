use FileUtil;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

generateClusters.pl - Take a set of RoC-ITS reads and cluster them into roughly genus clusters

=head1 VERSION

This document refers to version 1.00 of generateClusters.pl, released 07.06.2022.

=head1 SYNOPSIS

generateClusters.pl [-help]
Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message
    -reset                Overwrite existing results
    -noN                  Do not replace bases with N's during correction
    -input                Input fasta file of RoC-ITS sequences
    -outDir               Write final clusters to this directory
    -rocitsPath           Full path to where RoC-ITS scripts are located
                          default = "."
    -clusterPid           Clustering identity
                          default = 0.95 (designed but not guaranteed to cluster by genera)
    -minsubreads          Minimum number of subreads needed to process RoC-ITS sequence
    -minclusterSize       Minimum number of RoC-ITS sequences to make a valid genera cluster
                          default = 10
    -basename             Unique name that will identify the intermediate and final products of the processing
    -16s_fiveprime        Regular expression that identifies 5' end of 16S gene. Note that this is the reverse
                          complement of the 5' sequence.
                          default = "CTGAGCCAGGATCAAACTCT"
    -16s_threeprime       Regular expression that identifies 3' end of 16S gene. Note that this is the reverse
                          complement of the 3' sequence.
                          default = "TACGG.TACCTTGTTACGACTT"

=head1 DESCRIPTION

=head2 Overview

generateClusters.pl is designed to broadly automate the generation of genera-specific clusters from RoC-ITS reads

=cut

my ($help, $log, $msa, $output, $minCov, $excluded, $maxEdits, $reference, $noNs);
my ($rocitsPath);

BEGIN: {
  GetOptions(
    "help|h"            => \$help,
    "log"               => \$log,
    "input=s"           => \$input,
    "outDir=s"          => \$outdir,
    "basename=s"        => \$basename,
    "rocitsPath=s"      => \$rocitsPath,
    "minsubreads=s"     => \$minreads,
    "minclstrSize=s"    => \$minClusterSize,
    "16s_fiveprime=s"   => \$prime5,
    "clusterPid"        => \$cPid,
    "16s_threeprime=s"  => \$prime3,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}
if (!defined($input) || !defined($outdir) || !defined($basename) || !defined($minreads)) {
  pod2usage(1);
  exit;
}

$rocitsPath = "." unless defined($rocitsPath);
$prime3 = "TACGG.TACCTTGTTACGACTT" unless defined($prime3);
$prime5 = "CTGAGCCAGGATCAAACTCT" unless defined($prime5);
$minClusterSize = 10 unless defined($minClusterSize);
$cPid = 0.95 unless defined($cPid);

my $cmd;
$cmd = "perl $rocitsPath/selectFastaBasedOnSubreadCount.pl $input $minreads $rocitsPath > $basename.selected.fa";
print STDERR $cmd,"\n";
system($cmd);

$cmd = "perl $rocitsPath/identify16SBasedOnPrimerSequences.pl $basename.selected.fa $prime3 $prime5 $rocitsPath | perl $rocitsPath/reformatFasta.pl - 0 999999 > $basename.16Sonly.fa";
print STDERR $cmd,"\n";
system($cmd);

$cmd = "formatdb -p F -i $basename.16Sonly.fa";
print STDERR $cmd,"\n";
system($cmd);

$cmd = "blastall -i $input -d $basename.16Sonly.fa -p blastn -F F -a 20 -X 150 -q -5 -r 4 -v 5 -b 5 -e 1e-120 -m 8 > $basename.blastn.parsed";
print STDERR $cmd,"\n";
system($cmd);

$cmd = "perl $rocitsPath/getFullLength16SFromRoCITS.pl $basename.blastn.parsed $basename.16Sonly.fa $input $rocitsPath | sort > extract16S.$basename.sh";
print STDERR $cmd,"\n";
system($cmd);

$cmd = "bash extract16S.$basename.sh > $basename.16S.fa";
$cmd = "perl $rocitsPath/multiExtract.pl extract16S.$basename.sh $rocitsPath > $basename.16S.fa";
print STDERR $cmd,"\n";
system($cmd);

$cmd = "cd-hit-est -i $basename.16S.fa -o $basename.clustered.out -T 20 -M 0 -n 12 -c $cPid -d 5000 -l 20 -g 1";
print STDERR $cmd,"\n";
system($cmd);

$cmd = "mkdir -p $outdir";
print STDERR $cmd,"\n";
system($cmd);

$cmd = "perl $rocitsPath/processCDHitCluster.pl $basename.clustered.out.clstr $minClusterSize $rocitsPath | grep nt | perl $rocitsPath/produceClusters.pl - $basename.clustered.out.clstr $input $outdir/$basename $rocitsPath";
print STDERR $cmd,"\n";
system($cmd);
