use FileUtil;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Copyright (c) 2021, Indiana University, Bloomington. All rights reserved.

=head1 NAME

makeConsensus.pl - takes an MSA and identifies poorly represented regions on the ends of the MSA that
can be trimmed.

=head1 VERSION

This document refers to version 1.00 of processRoCITSSampleSet.pl, released 12.27.2021.

=head1 SYNOPSIS

makeConsensus.pl [-help]

Options:
    
    -log                  Print logging info to standard error
    -help|h               Print help message

    -msa                  Input fasta format MSA file
    -output               Output consensus fasta
    -rocitsPath           Full path to where RoC-ITS scripts are located

=head1 DESCRIPTION

=head2 Overview

Generates a consensus sequence from a MSA.

=cut

my ($help, $log, $fasta, $per);

BEGIN: {
  GetOptions(
    "help|h"        => \$help,
    "log"           => \$log,
    "msa=s"         => \$fasta,
    "output=s"      => \$output,
    "rocitsPath=s"  => \$rocitsPath,
  ) || pod2usage(2);
  pod2usage(1) if defined($help);
}

my $name = $fasta;
$name =~ s/.*\///;
$name =~ s/.fasta_aln$//;
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
		$totalDash += $dash;
		my @s = split //;
		my $l = @s;
		$len = $l if $l > $len;
		push @seq,\@s;
		push @ids,$id;
	}
}
$ffh->close();

my @con;
for (my $j = 0; $j < $len; $j++) {
	my %b;
	for (my $i = 0; $i < @ids; $i++) {
		$b{${$seq[$i]}[$j]}++;
	}
	foreach my $i (sort { $b{$b} <=> $b{$a};} keys %b) {
		push @con,$i;
		last;
	}
}
my $con = join "",@con;
$con =~ s/-//g;
$con =~ s/ //g;
my $perDash = int($totalDash/$totalChar*1000)/10;
my $fho = new FileHandle(">$output");
print $fho ">$name /count=$cnt /dash=$perDash\n";
print $fho $con,"\n";
$fho->close();
