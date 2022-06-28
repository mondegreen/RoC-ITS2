use FileUtil;

my $pick = shift;
my $cls = shift;
my $fasta = shift;
my $base = shift;
my $noCnt = shift;
my $centroidOnly = shift;

my %clstrs;
my $pfh = FileUtil::openFileHandle($pick);
while (<$pfh>) {
  chomp;
  if (/^(\d+)/) {
    $clstrs{$1} = 1;
  }
}
$pfh->close();

my $cfh = FileUtil::openFileHandle($cls);
while (<$cfh>) {
  if (/^>Cluster\s+(\d+)/) {
    $c = $1;
  } else {
    /^\S+\s+\S+\s+>(\S+)\.\.\./;
    my $id = $1;
    if ($centroidOnly) {
      if (/\.\.\. \*/) {
        $m{$c}{$id} = 1;
      }
    } else {
      $m{$c}{$id} = 1;
    }
  }
}
$cfh->close();

my $ffh = FileUtil::openFileHandle($fasta);
while (<$ffh>) {
  if (/^>(\S+)/) {
    my $id = $1;
    $_ = <$ffh>;
    chomp;
    $seq{$id} = $_;
  }
}
$ffh->close();

foreach my $i (keys %clstrs) {
  my $cnt = keys %{$m{$i}};
  my $fout;
  if (defined($noCnt) && $noCnt) {
    $fout = new FileHandle(">$base.$i.fa");
  } else {
    $fout = new FileHandle(">$base.$i.num=$cnt.fa");
  }
  foreach my $j (keys %{$m{$i}}) {
    print $fout ">$j\n$seq{$j}\n";
  }
  $fout->close();
}
