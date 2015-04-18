#!/usr/bin/perl
use strict;
open(INSPECT, $ARGV[0]) || die "cannot open inspect result file\n";

while(<INSPECT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\t/, $line);
    my $peptide = $tokens[2];
    $peptide =~ /^.\.(.+)\..$/;
    $peptide = $1;
    #print "peptide is $peptide\n";
    $peptide =~ /([\+\-]\d+)/;
    my $ptm = $1;
    if($ptm =~ /[0-9]/){
	print $ptm."\n";
    }
}
