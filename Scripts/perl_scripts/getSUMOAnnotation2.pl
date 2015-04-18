#!/usr/bin/perl
#print result file from MXDB crossslinke search to a annotation format

use strict;
open(RESULT, $ARGV[0]) || "die cannot open result file";

while(<RESULT>){
    my $line = $_;
    my @tokens = split(/\s+/, $line);
    $tokens[8] =~ s/[0-9\.\+]//g;
    if($tokens[5] =~ /\+\-17/){
	print "Scan Number: $tokens[3]\t".($tokens[8])."--".(substr($tokens[10], 0, 14))."\n";
    }else{
	print "Scan Number: $tokens[3]\t".(substr($tokens[5], 6, 12))."--".(substr($tokens[5], 0, 6))."\n";
    }
}
