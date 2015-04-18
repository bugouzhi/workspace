#!/usr/bin/perl
use strict;
open(MXDB, $ARGV[0]) || die "cannot open output file";
open(OUT, ">tempmxdb.txt");
my $forwarddb = $ARGV[1];
my $decoydb = $ARGV[2];
my $outfile = $ARGV[3];

while(<MXDB>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    $tokens[10] =~ s/[0-9\.\+]//g;
    $tokens[12] =~ s/[0-9\.\+]//g;
    for(my $i = 0; $i <= $#tokens; $i++){
	print OUT $tokens[$i]."\t";
    }
    print OUT "\n";
}
close(OUT);
`./mapDecoyMixturePeptides2.pl $forwarddb $decoydb tempmxdb.txt > $outfile`;
#`rm tempmxdb.txt`;

