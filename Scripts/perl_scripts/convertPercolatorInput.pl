#!/usr/bin/perl
#convert input file from mxdb to percolator format

use strict;
open(RESULT, $ARGV[0]) || die "cannot open input file";
my $id = 0;
print "PSMId\tLabel\tScannum\tScore\tScore1\tScore2\texplInt\tb1\ty1\tb2\ty2\tbseries1\tyseries1\tbseries2\tyseries2\tmerror1\tmerror2\texpInt1\texpInt2\tpMassDiff\tPeptide\tProteins\tScore\t-\t-\n";
while(<RESULT>){
    chomp $_;
    my @tokens = split(/\t/, $_);
    my $label = 1;
    if($tokens[5] =~ /X_/){
	$label = -1;
    }
    print "$id\t$label\t$tokens[1]\t";
    for(my $i = 7; $i < 25; $i++){
	if($i >= 7 & $i <= 8){
	    next;
	}
	print "$tokens[$i]\t";
    }
    my $precursorDiff = $tokens[2] - $tokens[6];
    if($precursorDiff < 0){
	$precursorDiff = -1*$precursorDiff;
    }
    #$precursorDiff = $precursorDiff*1000000/$tokens[2];
    print "$precursorDiff\t";
    $id++;
    my $pep = $tokens[4];
    my @tokens2 = split(" & ", $pep);
    my $pep = "K.$tokens2[0].K";
    my $prot = $tokens[5];
    #@tokens2 = split(/ & /, $prot);
    $prot =~ s/\s/_/g;
    #print
    print "$pep\t$prot\n";
}
