#!/usr/bin/perl
use strict;
open(RESULT, $ARGV[0]);

while(<RESULT>){
    my $line = $_;
    chomp $line;
    if($line =~ /best/){
	my @tokens = split(/\s+/, $line);
	for(my $i = 0; $i < 6; $i++){
	    print $tokens[$i]."\t";
	}
	for(my $i = 6; $i < $#tokens; $i++){
	    printf("%.2f\t", $tokens[$i]);
	}
	printf("%.2f\n", $tokens[$#tokens]);
    }
}
