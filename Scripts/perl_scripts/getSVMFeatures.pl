#!/usr/bin/perl
#extract SVM features from MXDB results file
use strict;

open(RESULT, $ARGV[0]) || die "cannot open result file\n";
my $beginIndex = $ARGV[1];
my $label = $ARGV[2];

while(<RESULT>){
    my $line = $_;
    #print "line is $line\n";
    chomp $line;
    my @tokens = split(/\s+/, $line);
    print "$label ";
    for(my $i = $beginIndex; $i < $#tokens; $i++){
	print ($i-$beginIndex+1);
	print":".$tokens[$i]."\t";
    }
    print $#tokens-$beginIndex+1;
    print ":".$tokens[$#tokens]."\n";
}


