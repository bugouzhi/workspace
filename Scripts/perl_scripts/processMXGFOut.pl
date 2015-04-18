#!/usr/bin/perl
use strict;

open(RESULT, $ARGV[0]) || die "cannot open result files";

while(<RESULT>){
    my $line = $_;
    chomp $line;
    my $line2 = <RESULT>;
    chomp $line2;
    my @tokens = split(/\t/, $line);
    my @tokens2 = split(/\t/, $line2);
    my $probaInd = 36;
    my $ind1 = 36;
    my $ind2 = 37;
    my $prot1 = 7;
    my $prot2 = 8;
    my $pep1 = 5;
    my $pep2 = 6;
    if($tokens2[$probaInd] >= $tokens[$probaInd]){#&& ($tokens2[35] < $tokens2[36] || $tokens[34] > $tokens2[34])){
	for(my $i = 0; $i <= $#tokens; $i++){
	    print $tokens[$i]."\t";
	}
	print "\n";
    }else{
	$tokens2[$ind1] = $tokens[$ind2];
	$tokens2[$ind2] = $tokens[$ind1];
	$tokens2[$prot1] = $tokens[$prot2];
	$tokens2[$prot2] = $tokens[$prot1];
	$tokens2[$pep1] = $tokens[$pep2];
	$tokens2[$pep2] = $tokens[$pep1];
	for(my $i = 0; $i <= $#tokens2; $i++){
	    print $tokens2[$i]."\t";
	}
	print "\n";
    }
    if($tokens[3] =~ /^r/){
	#$tokens2[2] = "r".$tokens2[2];
    }
    if($tokens2[3] =~ /^r/){
	#$tokens[2] = "r".$tokens[2];
    }
    if($tokens[4] > $tokens[5]){
	for(my $i = 0; $i < $#tokens; $i++){
	    #print $tokens[$i]."\t";
	}
	#print "\n";

    }else{
	if($tokens2[4] > $tokens[5]){
	    for(my $i = 0; $i < $#tokens2; $i++){
		#print $tokens2[$i]."\t";
	    }
	    print "\n";
	}else{
	    for(my $i = 0; $i < $#tokens; $i++){
		#print $tokens[$i]."\t";
	    }
	    #print "\n";
	}
    }
}
