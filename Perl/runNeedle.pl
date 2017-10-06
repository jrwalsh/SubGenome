use strict;
use warnings;

my ($aseq, $bseq) = @ARGV;

my @args = ("bash", "-c", "needle <(echo $aseq) <(echo $bseq) -gapopen 10 -gapextend 0.5 -outfile >(cat)");
system(@args);
