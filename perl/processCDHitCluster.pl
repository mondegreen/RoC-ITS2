use FileUtil;

my $f = shift;
my $mmin = shift;
my $per_in = shift;
my $keyword = shift;

my $fh = FileUtil::openFileHandle($f);
while (<$fh>) {
  if (/^>Cluster\s+(\d+)/) {
    $c = $1;
  } else {
    $cnt{$c}++;
    $total++;
  }
}
$fh->close();

my $min = 99999999;
my $fh = FileUtil::openFileHandle($f);
while (<$fh>) {
  if (/^>Cluster\s+(\d+)/) {
    $c = $1;
  } else {
    if (/$keyword/o) {
      $min = $cnt{$c} if $cnt{$c} < $min;
    }
  }
}
$fh->close();

if ($mmin && !defined($keyword)) {
  $min = $mmin;
}
print $min,"\n";
my $fh = FileUtil::openFileHandle($f);
while (<$fh>) {
  if (/^>Cluster\s+(\d+)/) {
    $c = $1;
    if ($cnt{$c} >= $min) {
      $seen += $cnt{$c};
    }
  } else {
    if ($cnt{$c} >= $min) {
      print "$c\t$cnt{$c}\t$_";
    }
  }
}
$fh->close();

my $per = int(1000*$seen/$total)/10;
print "Incorp $per_in $total $seen $per\n";
