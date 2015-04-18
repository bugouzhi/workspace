#!/usr/bin/perl
use strict;
open(LIST1, $ARGV[0]) || die "cannot open file1";
open(LIST2, $ARGV[1]) || die "cannot open file2";

my %pepTable = ();
while(<LIST1>){
	my $line = $_;
	chomp $line;
	my @tokens = split(/\t/, $line);
	my @tokens2 = split(/\./, $tokens[2]);
	#print "storing $tokens2[1]\n";
	$pepTable{$tokens2[1]} = $tokens[1];
}

while(<LIST2>){
	my $line = $_;
	chomp $line;
	my @tokens = split(/\./, $line);
	if(exists $pepTable{$tokens[0]}){
		print "found peptide $line in both\n";
	}
}
