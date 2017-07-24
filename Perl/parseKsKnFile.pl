use strict;
use warnings;

my $block = "";
print "block\tks\tkn\torg_chr1\tstart1\tstop1\tgene1\torg_chr2\tstart2\tstop2\tgene2\n";
while(<>) {
  my @words = split /\t/, $_;
  if ($words[0] =~ /#\d+/) {
    if ($words[4] =~ /f/) {
      $block = $words[0] . "f";
    } else {
      $block = $words[0] . "r";
    }
    #$block = $words[0];
    #print "$block \n";
  } elsif ($words[0] =~ /#Ks/ | $words[0] =~ /#This.*/) {
    #skip
  } else {
    my $ks = $words[0];
    my $kn = $words[1];
    my $org_chr1 = $words[2];
    my @org_features1 = split(/\|\|/, $words[3]);
    my $gene1 = $org_features1[3];
    my $start1 = $words[4];
    my $stop1 = $words[5];
    my $org_chr2 = $words[6];
    my @org_features2 = split(/\|\|/, $words[7]);
    my $gene2 = $org_features2[3];
    my $start2 = $words[8];
    my $stop2 = $words[9];
    my $eval = $words[10];
    my $block_score = $words[11];
    my $gevo_link = $words[12];

    print $block, "\t", $ks, "\t", $kn, "\t", $org_chr1, "\t", $start1, "\t", $stop1, "\t", $gene1, "\t", $org_chr2, "\t", $start2, "\t", $stop2, "\t", $gene2, "\n"; 
  }
}
