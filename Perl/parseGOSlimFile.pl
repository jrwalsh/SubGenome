use strict;
use warnings;

while(<>) {
  if ($_ =~ /id: GO:.*$/) {
    my @words = split / /, $_;
    print $words[1];
  }
}
