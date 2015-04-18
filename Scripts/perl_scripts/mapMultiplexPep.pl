#!/usr/bin/perl
use strict;
open(PEPS, $ARGV[0]) || die "cannot open peptides file\n";
open(RESULT, $ARGV[1]) || die "cannot open result file";

my %pepTable={};
while(<PEPS>){
    my $pep = $_;
    chomp $pep;
    $pepTable{$pep} = 1;
}

while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\t/, $line);
    my @peps = split(/ & /, $tokens[5]);
    for(my $i = 0; $i < 3; $i++){
	my $pep = $peps[$i];
	$pep = substr($pep, 0, length($pep)-2);
	#print "peptide is $pep\n";
	if(!$pepTable{$pep}){
	    print "$tokens[1]\t$pep not found\n"
	}
    }
}
