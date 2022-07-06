use FileUtil;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

processRoCITSSampleSet.pl - Takes a set of genus specific reads, cluster them at a specific percent
  identity cutoff, and generate an MSA and consensus sequence for each cluster.

=head1 VERSION

This document refers to version 1.00 of processRoCITSSampleSet.pl, released 12.27.2021.

=head1 SYNOPSIS

processRoCITSSampleSet.pl [-help]

Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message
    -reset                Overwrite existing results

    -fasta                Input fasta file
    -targerDir            Directory to create and store results in
    -minSize              Minimum size of clusters to keep
    -overrideMinCnt       
    -pid                  Clustering percent identity (times 10 so 970 equals 97.0%)
    -identifier           Identifier to uniquely describe sample intermediate files
    -jobs                 Max number of concurrent jobs
    -chimeraDb            File containing the trusted reference sequences
    -rocitsPath           Full path to where RoC-ITS scripts are located

    -reference            (Optional) Reference indicating safe identifiers
    -muscle               (Optional) Full path to Muscle executable
    -probcons             (Optional) Full path to probcons executable
    -vsearch              (Optional) Full path to vsearch executable
    -context              (Optional) Region around variable bases to keep (default 2)

=head1 DESCRIPTION

=head2 Overview

processRoCITSSampleSet.pl executes a pipeline to cluster RoC-ITS reads at a specified
percent identity.

=cut

my ($help, $log, $reset, $input, $workDir, $minClstPercent, $minClst, $overrideMinCnt, $reference, $chimeradb, $id, $jobs);
my ($muscleBin, $probconsBin, $rocitsPath, $context, $vsearchBin);

BEGIN: {
  GetOptions(
    "help|h"            => \$help,
    "log"               => \$log,
    "reset"             => \$reset,
    "fasta=s"           => \$input,
    "targetDir=s"       => \$workDir,
    "pid=s"             => \$minClstPercent,
    "minSize=s"         => \$minClst,
    "context=s"         => \$context,
    "overrideMinCnt:s"  => \$overrideMinCnt,
    "chimeraDb=s"       => \$chimeradb,
    "identifier=s"      => \$id,
    "maxJobs:s"         => \$jobs,
    "rocitsPath=s"      => \$rocitsPath,
    "reference:s"       => \$reference,
    "muscle:s"          => \$muscleBin,
    "vsearch:s"         => \$vsearchBin,
    "probcons:s"        => \$probconsBin,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}

$reference = "NotAnIdentifierX" unless defined($reference);
$jobs = 40 unless defined($jobs);
$muscleBin = "muscle" unless defined($muscleBin);
$vsearchBin = "vsearch" unless defined($vsearchBin);
$probconsBin = "probcons" unless defined($probconsBin);
$minClstPer = $minClstPercent;
while ($minClstPer > 1) {
  $minClstPer = $minClstPer/10;
}
$maxClstPer = $maxClstPercent;
while ($maxClstPer > 1) {
  $maxClstPer = $maxClstPer/10;
}
$context = 2 unless defined($context);

my $cmd;
my $loc = `pwd`;
chomp $loc;

my $id2 = "$id.$minClstPercent";

$cmd = "mkdir -p $workDir";
runCommand($cmd, $log, __LINE__);

### Perform initial multiple sequence alignment w/ muscle
my $vsearchOut = "$workDir/$id.noChim.fa";
if ($chimeradb) {
  $cmd = "$vsearchBin --db $chimeradb --uchime_ref $input --nonchimeras $vsearchOut";
} else {
  $cmd = "cp $input $vsearchOut";
}
my $size = -s $vsearchOut;
if ($reset || !$size) {
  runCommand($cmd, $log, __LINE__);
}

### Perform initial multiple sequence alignment w/ muscle
my $muscleOut = "$workDir/$id.musc.fa";
$cmd = "$muscleBin -super5 $vsearchOut -threads $jobs -output $muscleOut";
my $size = -s $muscleOut;
if ($reset || !$size) {
  runCommand($cmd, $log, __LINE__);
}

### determine number of input sequences
$num = `grep -c \"^>\" $input`;
chomp $num;
if ($log) {
  print STDERR "LOG: Input reads = $num: $0",__LINE__,"\n";
}

### muscle MSA's have certain common errors that can be corrected with the improveMSA.pl script
### Calculate the minimum coverage to keep a column in MSA. The value will be either 3 or 0.
### The value 3 is used if the number of sequences * 5% is greater than 3.
my $minCnt = int($num * 0.05);
$minCnt = $minCnt > 3 ? 3:0;
my $improvedOut = "$workDir/$id.musc.improved.fa";
$cmd = "perl $rocitsPath/improveMSA.pl -log $log -maxJobs $jobs -msa $muscleOut -wordSize 12 -minimumCov $minCnt -probcons $probconsBin -local $workDir -rocitsPath $rocitsPath -out $improvedOut";
$size = -s $improvedOut;
if ($reset || !$size) {
  runCommand($cmd, $log, __LINE__);
}

### correct the MSA to eliminate error/noise in the data. This is very valuable for isolates but
### can be problematic for complex mixed collections of organisms. 
my $minCnt = int($num * 0.05);
$minCnt = $minCnt < 3 ? 3:$minCnt;
$minCnt = $overrideMinCnt if defined($overrideMinCnt);
my $correctedOut = "$workDir/$id.musc.corrected.fa";
my $excludedOut = "$workDir/$id.musc.excluded.fa";
$cmd = "perl $rocitsPath/correctMSA.pl -rocitsPath $rocitsPath -log $log -msa $improvedOut -minCov $minCnt -reference $reference -out $correctedOut -excluded $excludedOut";
$size = -s $correctedOut;
if ($reset || !$size) {
  print "LOG: Minimum coverage for correction = $minCnt\n" if $log;
  runCommand($cmd, $log, __LINE__);
}

$cmd = "perl $rocitsPath/findBoundaries.pl -msa $correctedOut";
my $bounds = execCommand($cmd, $log, __LINE__);
chomp $bounds;
my @bounds = split /\s+/,$bounds;
$bounds[0] += 10;
$bounds[1] += 10;
print "LOG: Left bound = $bounds[0] Right bound = $bounds[1]: $0 line ",__LINE__,"\n";

my ($variantRegionsOut,$functionalRegionsOut,$cdhitOut);
if (defined($bounds)) {
  $variantRegionsOut = "$workDir/$id.interesting.$context.regions.fa";
  $cmd = "perl $rocitsPath/identifyVariantColumns.pl -rocitsPath $rocitsPath -log $log -msa $correctedOut -pid $minClstPer -leftTrim $bounds[0] -rightTrim $bounds[1] -context $context -minColumns 20 -output $variantRegionsOut";
  $size = -s $variantRegionsOut;
  if ($reset || !$size) {
    runCommand($cmd, $log, __LINE__);
  }

  $functionalRegionsOut = "$workDir/$id.functional.regions.fa";
  $cmd = "perl $rocitsPath/extractEffectiveColumns.pl -rocitsPath $rocitsPath -log $log -msa $correctedOut -leftTrim $bounds[0] -rightTrim $bounds[1] -output $functionalRegionsOut";
  $size = -s $functionalRegionsOut;
  if ($reset || !$size) {
    runCommand($cmd, $log, __LINE__);
  }

  $cdhitOut = "$workDir/$id2.cd-hit.cont.$context.regions.fa";
  $cmd = "cd-hit-est -i $variantRegionsOut -o $cdhitOut -T 20 -M 0 -n 5 -c $minClstPer -d 5000 -l 10 -g 1";
  $size = -s $cdhitOut;
  if ($reset || !$size) {
    runCommand($cmd, $log, __LINE__);
  }
} else {
  exit;
}

exit unless defined($minClst);

$targetDir = "$workDir/$id2.cont.$context.clusters";
if (-e "$targetDir") {
  $cmd = "rm -r $targetDir";
  runCommand($cmd, $log);
}
$cmd = "mkdir -p $targetDir";
runCommand($cmd, $log);
print STDERR "LOG: $cmd: $0 line ",__LINE__,"\n" if $log;
chdir "$targetDir";
print STDERR "LOG: chdir $workDir/$targetDir: $0 line ",__LINE__,"\n" if $log;


my $minCnt = int($num * 0.05);
$minCnt = $minClst if $minClst;
$cmd = "perl /N/projects/rollingITS/drusch/RocSplint/subgroup8ecoli/processCDHitCluster2.pl $loc/$id2.cd-hit.regions.fa.clstr $minCnt -reference $reference | grep nt | perl /N/projects/rollingITS/drusch/RocSplint/subgroup8ecoli/extractClusters.pl - $loc/$id.functional.regions.fa > $loc/$id.selected.functional.fa";
print STDERR $cmd,"\n";
#system($cmd);

$cmd = "perl $rocitsPath/associateSequencesWithClusters.pl -rocitsPath $rocitsPath -log $log -clstr $loc/$cdhitOut.clstr -minSize $minCnt | grep nt | perl $rocitsPath/createSubclusters.pl -rocitsPath $rocitsPath -log $log -selected - -clstr $loc/$cdhitOut.clstr -fasta $loc/$correctedOut -baseName picked.$id2.clusters -addCount 1";
runCommand($cmd, $log, __LINE__);

#$cmd = "perl /N/projects/rollingITS/drusch/RocSplint/subgroup8ecoli/processCDHitCluster2.pl $loc/$id2.cd-hit.regions.fa.clstr $minCnt | grep nt | perl /N/projects/rollingITS/drusch/RocSplint/subgroup8ecoli/selectClusters.pl - $loc/$id2.cd-hit.regions.fa.clstr $loc/$id.functional.regions.fa picked.clusters 1";
#print STDERR $cmd,"\n";
#system($cmd);
#exit;

opendir my $dir, "." or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;

my @ids;
foreach my $i (@files) {
  if ($i =~ /(picked.$id2.clusters.(\d+).fa)/) {
    push @ids,$2;
    print STDERR "LOG: Identified file $1 for processing: $0 line ",__LINE__,"\n" if $log;
  }
}
#exit;

my $list = join " ",sort {$a <=> $b;} @ids;
my $ofh = new FileHandle(">process.sh");
print $ofh "setenv PATH $ENV{PATH}\n";
print $ofh "foreach I ( $list )\n";
print $ofh "  $muscleBin -super5 picked.$id2.clusters.\$I.fa -threads $jobs -output picked.$id2.clusters.\$I.muscle.fa; perl $rocitsPath/improveMSA.pl -log $log -maxJobs $jobs -msa picked.$id2.clusters.\$I.muscle.fa -wordSize 8 -minimumCov 3 -probcons $probconsBin -local . -rocitsPath $rocitsPath -out picked.$id2.clusters.\$I.muscle.improved.fa;  ; perl $rocitsPath/makeConsensus.pl -msa picked.$id2.clusters.\$I.muscle.improved.fa -rocitsPath $rocitsPath -output picked.$id2.clusters.\$I.muscle.improved.consensus.fa\n";
print $ofh "end\n";
$ofh->close();
$cmd = "$muscleBin -super5 $vsearchOut -threads $jobs -output $muscleOut";

$cmd = "tcsh -x process.sh";
runCommand($cmd, $log, __LINE__);
#print STDERR $cmd,"\n";
#system($cmd);

exit;

sub runCommand {
  my ($cmd,$state,$pos) = @_;
  
  if ($state) {
    print STDERR "LOG: $cmd: $0 line $pos\n";
  }
  my $result = system($cmd);
  if ($result) {
    my $err = $? >> 8;
    print STDERR "ERROR: $err: $0 line ",__LINE__,"\n";
  }
}

sub execCommand {
  my ($cmd,$state,$pos) = @_;
  
  if ($state) {
    print STDERR "LOG: $cmd: $0 line $pos\n";
  }
  my $result = `$cmd`;
  return $result;
}
