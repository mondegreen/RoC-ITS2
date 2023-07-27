use FileUtil;

my $f = shift;
my $mol = shift;
my $b = shift;
my $e = shift;

my $fh = FileUtil::openFileHandle($f);
while (<$fh>) {
  if (/^>(\S+)/) {
    $newid = $1;
    if (@seq) {
      $seq{$id} = join "",@seq;
      undef @seq;
    }
    $id = $newid;
    next;
  } else {
    chomp;
    push @seq,$_;
  }
}
if (@seq) {
  $seq{$id} = join "",@seq;
}
my $sub = substr($seq{$mol},$b,$e-$b);
my $len = length($sub);
my $full = length($seq{$mol});
print ">$mol /length=$len /full_length=$full /begin=$b /end=$e\n$sub\n";
