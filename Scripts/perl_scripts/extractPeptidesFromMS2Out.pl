#!/usr/bin/perl
#given a list of peptide, remove any result that is not part of the search space for MDBSearches
use strict;
open(SEARCHSPACE, $ARGV[0]) || "die cannot open peptide list file";
open(RESULTFILE, $ARGV[1]) || "die cannot open search result file ";
my %table=();
while(<SEARCHSPACE>){
    my $line = $_;
    chomp $line;
    $table{$line} = 1;
}

while(<RESULTFILE>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    my @tokens2 = split(/\./, $tokens[3]);
    #print "id is $tokens[3]\n";
    if(exists $table{$tokens2[1]}){
	print $line."\n";
    }
}
