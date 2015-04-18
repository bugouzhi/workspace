#!/usr/bin/perl
#print result file from MXDB single peptide search to a annotation format

use strict;
open(RESULT, $ARGV[0]) || "die cannot open result file";

while(<RESULT>){
    my $line = $_;
    my @tokens = split(/\s+/, $line);
    if($tokens[5] =~ /Scan/){
	$tokens[10] =~ /(\+5\d\d\.[0-9\.\+]+)/;
	my $offset = $1;
	#print "offset is $offset\n";
	my $position = index($tokens[10], $offset);
	$tokens[10] =~ s/\+5\d\d\.[0-9\.\+]+//g;
	my $mod = index($tokens[10], "+"); #look for additional mods
	#print "mod is $mod\n";
	if($mod >= 0 && $mod < $position){
	    $position = $position - 7;
	}
	$tokens[12] =~ /(Q[0-9\.\-]*QQTGG)/;
	my $tag = $1;
		    print "Scan Number: $tokens[3]\t".($tokens[10])."--".$tag."\t".$position."\t6\t-18.0106\n";		    
    }else{
	if($tokens[5] =~ /(Q[0-9\-\.\+]*QQTGG)(K[0-9\+\.]*A.+)/){
	    print "Scan Number: $tokens[3]\t".$2."--".$1."\n";
	}
    }
}
