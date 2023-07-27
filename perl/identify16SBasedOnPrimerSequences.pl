use FileUtil;

my $rocITSFasta = shift;
my $primerA = shift;
my $primerB = shift;
my $rocitsPath = shift;

### default primers to identify 16S genes in eubacteria
$primerA = "TACGG.TACCTTGTTACGACTT" unless defined($primerA);
$primerB = "CTGAGCCAGGATCAAACTCT" unless defined($primerB);

my $rfh = FileUtil::openFileHandle("perl $rocitsPath/reformatFasta.pl $rocITSFasta 0 999999999|");
while (<$rfh>) {
  if (/^>/) {
    $def = $_;
  } else {
    if (/.*($primerA.*$primerB)/o) {
      print $def;
      print $1,"\n";
    }
  }
}
$rfh->close();
