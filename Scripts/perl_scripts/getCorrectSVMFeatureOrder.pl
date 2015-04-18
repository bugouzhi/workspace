#!/usr/bin/perl
#correct the svm feature orders for MXDB output

use strict;

open(OUT, $ARGV[0]) || die "cannot open database output file";

while(<OUT>){
    my $line = $_;
    if($line =~ /best:/){
	chomp $line;
	my @tokens = split(/\s+/, $line);
	if($tokens[15] >= $tokens[16]){
	    for(my $i = 0; $i < 32; $i++){
		print $tokens[$i]."\t";
	    }
	    print "\n";
	    #print $tokens[$#tokens]."\n";
	}else{
	    my $pep = $tokens[12];
	    $tokens[12] = $tokens[10];
	    $tokens[10] = $pep;
	    for(my $i = 0; $i < 15; $i++){
		print $tokens[$i]."\t";
	    }
	    print $tokens[16]."\t".$tokens[15]."\t".$tokens[18]."\t".$tokens[17]."\t";
	    print $tokens[19]."\t";
	    print $tokens[22]."\t".$tokens[23]."\t".$tokens[20]."\t".$tokens[21]."\t";
	    print $tokens[26]."\t".$tokens[27]."\t".$tokens[24]."\t".$tokens[25]."\t";
	    print $tokens[29]."\t".$tokens[28]."\t".$tokens[31]."\t".$tokens[30]."\n";
	}
    }
}

