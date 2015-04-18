#!/usr/bin/perl
use strict;
#generate all possible peptide sequene from a protein sequence

open(FASTA, $ARGV[0]) || die "cannot open fasta files\n";
my $minLength = $ARGV[1];
my $maxLength = $ARGV[2];

while(<FASTA>){
    if($_ =~ />/){
	#print "line is $_\n";
	my $proteinSeq = <FASTA>;
	chomp $proteinSeq;
	for(my $begin = 0; $begin < length($proteinSeq); $begin++){
	    for(my $length = $minLength; $length <= $maxLength; $length++){
		if($begin+$length <= length($proteinSeq)){
		    #print $begin."\t".$length."\n";
		    print substr($proteinSeq, $begin, $length)."\n";;
		}
	    }
	}
    }
}
