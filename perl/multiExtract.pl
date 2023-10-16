use FileUtil;

my $extract = shift;
my $rocitsPath = shift;

my %files;
my $efh = FileUtil::openFileHandle($extract);
while (<$efh>) {
  chomp;
  my @d = split /\s+/;
  $files{$d[2]} = 1;
}
$efh->close();

foreach my $i (keys %files) {
  my $fh = FileUtil::openFileHandle("perl $rocitsPath/reformatFasta.pl $i 0 999999|");
  while (<$fh>) {
    chomp;
    if (/^>(\S+)/) {
      my $id = $1;
      if (exists($def{$id})) {
        print STDERR "Duplicate $id\n";
        next;
      }
      $def{$id} = $_;
      $_ = <$fh>;
      chomp;
      $seq{$id} = $_;
      $len{$id} = length($_);
    }
  }
  $fh->close();
}

my $efh = FileUtil::openFileHandle($extract);
while (<$efh>) {
  chomp;
  my @d = split /\s+/;
  my $sub = substr($seq{$d[3]},$d[4],$d[5]-$d[4]);
  my $len = length($sub);
  print ">$d[3] /length=$len /full_length=$len{$d[3]} /begin=$d[4] /end=$d[5]\n$sub\n";
}
$efh->close();
