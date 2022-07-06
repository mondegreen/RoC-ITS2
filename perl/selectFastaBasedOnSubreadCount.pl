use FileUtil;

my $fasta = shift;
my $cnt = shift;

my $fh = FileUtil::openFileHandle($fasta);
while (<$fh>) {
  if (/^>/) {
    if (/subreads=(\d+)/) {
      if ($1 >= $cnt) {
        print $_;
        $_ = <$fh>;
        print $_;
      }
    }
  }
}
$fh->close();
