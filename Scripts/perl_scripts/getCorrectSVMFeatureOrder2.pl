#!/usr/bin/perl
#correct the svm feature orders for MXDB output

use strict;

open(OUT, $ARGV[0]) || die "cannot open database output file";

while(<OUT>){
    my $line = $_;
    if($line =~ /best/){
	chomp $line;
	my @tokens = split(/\s+/, $line);
	if($tokens[14] >= $tokens[15]){
	    for(my $i = 0; $i < $#tokens; $i++){
		print $tokens[$i]."\t";
	    }
	    print $tokens[$#tokens]."\n";
	}else{
	    for(my $i = 0; $i < 14; $i++){
		print $tokens[$i]."\t";
	    }
	    print $tokens[15]."\t".$tokens[14]."\t".$tokens[17]."\t".$tokens[16]."\t";
	    print $tokens[18]."\t";
	    print $tokens[21]."\t".$tokens[22]."\t".$tokens[19]."\t".$tokens[20]."\t";
	    print $tokens[25]."\t".$tokens[26]."\t".$tokens[23]."\t".$tokens[24]."\t";
	    print $tokens[28]."\t".$tokens[27]."\t".$tokens[30]."\t".$tokens[29]."\t".$tokens[31]."\t".$tokens[32]."\n";
	}
    }
}

