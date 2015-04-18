#!/usr/bin/perl
#convert inspect output to linked peptides 
use strict;
open(INSPECT, $ARGV[0]) || die "cannot open output file";
open(OUT, ">temp_mapInspectToLinked.txt");

while(<INSPECT>){
    my $line = $_;
    chomp $line; 
    my @tokens = split(/\t/, $line);
    print "Scan\tNumber:\t$tokens[1]\t".substr($tokens[2],8,12)."--".substr($tokens[2],2,6)."\t".$tokens[22]."\n";
}
