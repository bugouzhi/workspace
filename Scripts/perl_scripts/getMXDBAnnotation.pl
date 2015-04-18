#!/usr/bin/perl
#print result file from MXDB single peptide search to a annotation format

use strict;
open(RESULT, $ARGV[0]) || "die cannot open result file";

my $scan=1;
my $peptide1=10;
my $peptide2=12;
my $offSetMass=$ARGV[1];
while(<RESULT>){
    my $line = $_;
    my @tokens = split(/\s+/, $line);
    if($tokens[$scan] =~ /Scan/){
	$tokens[$peptide1] =~ /([\+\-]\d\d\d+\.[0-9\.\+]+)/;
	my $offset = $1;
	#print "offset is $offset\n";
	my $position1 = index($tokens[$peptide1], $offset);
	my $subpep1 = substr($tokens[$peptide1], 0, $position1);
	$subpep1 =~ s/[0-9\.\+\-]//g;
	$position1 = length($subpep1);
	$tokens[$peptide1] =~ s/([\+\-]\d\d\d+\.[0-9\.\+\.]+)//g;
	$tokens[$peptide2] =~ /([\+\-]\d\d\d+\.[0-9\.\+\.]+)/;
	my $offset2 = $1;
	#print "offset is $offset\n";
	my $position2 = index($tokens[$peptide2], $offset2);
	my $subpep2 = substr($tokens[$peptide2], 0, $position2);
	$subpep2 =~ s/[0-9\.\+\-]//g;
	$position2 = length($subpep2);
	$tokens[$peptide2] =~ s/([\+\-]\d\d\d+\.[0-9\.\+\.]+)//g;
	print "Scan Number: $tokens[3]\t".($tokens[$peptide1])."--".($tokens[$peptide2])."\t".$position1."\t".$position2."\t$offSetMass\n";		    
	#print "Scan Number: $tokens[3]\t".($tokens[$peptide1])."--".($tokens[$peptide2])."\t7\t7\t138.069\n";		    
    }else{
	if($tokens[5] =~ /(Q[0-9\-\.\+]*QQTGG)(K[0-9\+\.]*A.+)/){
	    print "Scan Number: $tokens[3]\t".$2."--".$1."\n";
	}
    }
}
