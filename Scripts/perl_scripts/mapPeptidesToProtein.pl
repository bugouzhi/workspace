#!/usr/bin/perl
use strict;
open(PROT, $ARGV[0]) || die "cannot open protein file";
open(PEP, $ARGV[1]) || die "cannot open peptide file";
my $index = $ARGV[2];

my %table = ();

my $header = "";
my $currentProt="";
while(<PROT>){
    my $line = $_;
    chomp $line;
    if($line =~ />/){
	$table{$header} = $currentProt;
	$header = $line;
	$currentProt = "";
    }else{
	$currentProt = $currentProt.$line;
    }
}
#print "Done loading proteins\n";


$currentProt = "";
while(<PEP>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    foreach my $key(keys %table){
	#print "key is $key\n";
	$currentProt = $table{$key};
	#print "this prot: ".$currentProt."\n";
	#print "pep is $tokens[$index]\n";
	if($currentProt =~ $tokens[$index]){
	    print $line."\t".$key."\n";
	}
    }
}
