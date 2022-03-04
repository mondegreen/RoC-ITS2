package FileUtil;

use FileHandle;
use Carp;

sub openFileHandle {
  my ($file) = @_;
  $file = "-" unless defined($file);

  my $fh = new FileHandle();
  if ($file =~ /gz$/) {
     $fh->open("gunzip -c $file |") ||
      confess("ERROR: Problem opening $file");
  } elsif ($file =~ /bz2$/) {
     $fh->open("bunzip2 -c $file |") ||
      confess("ERROR: Problem opening $file");
  } elsif ($file eq "-") {
    $fh = $fh->fdopen(fileno(STDIN),"r") ||
      confess("ERROR: Problem opening STDIN");
  } else {
     $fh->open("$file") ||
      confess("ERROR: Problem opening $file");
  }

  return $fh;
}

sub getFileByFTP {
  my ($URL, $filename) = @_;
  
  if (!defined($filename)) {
    ($filename) = $URL =~ /.*\/(.*)/;
  }

use LWP::UserAgent;
  my $ua = LWP::UserAgent->new;

  my $req = HTTP::Request->new(GET => $URL);
  my $res = $ua->request($req, $filename);
  if ($res->is_success) {
  } else {
    # print any error messages
    print $res->status_line, "\n";
  }
}

1;
