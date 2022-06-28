use FileUtil;

my $tabBlastFile = shift;
my $reference16SFile = shift;
my $source16SFile = shift;

my $tfh = FileUtil::openFileHandle($tabBlastFile);
my $rfh = FileUtil::openFileHandle("perl reformatFasta.pl $reference16SFile |");
while (<$rfh>) {
  if (/^>(\S+).*full_length=(\d+)/) {
    $lengths{$1} = $2;
  }
}
$rfh->close();

while (<$tfh>) {
  chomp;
  my @d = split /\t/;
  if (exists($m{$d[0]})) {
  } else {
    $m{$d[0]} = $d[6]-1; ### -1 it to shift from base to inter-base coordinates
    $mm{$d[0]} = $d[7];
    $mmm{$d[0]} = $lengths{$d[1]}; ### associate length of the 16S sequence with the best alignment
  }
}

foreach my $i (keys %m) {
  my $l = join " ","perl getFastaSubset.pl",$source16SFile,$i,$m{$i},$mm{$i}; $x = ($mm{$i}-$m{$i})/$mmm{$i}; print $l,"\n" if $x >= 0.98;
}
