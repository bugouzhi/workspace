#!/usr/bin/perl
use strict;
#digest protein in silico into peptides

open(FASTA, $ARGV[0]) || die "cannot open fasta protein file\n";
my $maxLength = 45;
my $minLength = 7;
my $proteins = "*";
while(<FASTA>){
    my $line = $_;
    chomp $line;
    $line =~ s/\s+//g;
    if($line =~ />/){
	$proteins = $proteins."*";
    }else{
	$proteins = $proteins.$line;
    }
}

for(my $i = 0; $i < length($proteins); $i++){
    for(my $j = $minLength+2; $j <= $maxLength+2; $j++){
	my $peptide = substr($proteins, $i, $j+1);
	#if($peptide =~ /^[FYWLM]([A-Z]+[FYWLM])$/){
	#if($peptide =~ /^[FYWL*]([A-Z]+[FYWL][A-Z])$/ || $peptide =~ /^[FYWL][A-Z]+\*$/){
	#if($peptide =~ /^[A-Z]([A-Z]+[RK])$/ || $peptide=~ /^[KR][A-Z]+$/){
	if($peptide =~ /^[KR*]([A-Z]+[RK][A-Z])$/ || $peptide =~ /^[KR][A-Z]+\*$/){
	#if($peptide =~ /^[EQ*]([A-Z]+[EQ][A-Z])$/ || $peptide =~ /^[EQ][A-Z]+\*$/){
	#if($peptide =~ /^[A-Z][DEKR][A-Z]+[DEKR*]$/ || $peptide =~ /^\*[A-Z]+[DEKR*]$/){
	#if($peptide =~ /^[KR*][A-Z]+/ || $peptide =~ /^[A-Z]+[KR][A-Z]$/){
	#if($peptide =~ /^[A-Z*]([A-Z]+[A-Z][A-Z])$/ || $peptide =~ /^[A-Z][A-Z]+\*$/){
	#if($peptide =~ /^[A-Z]([ED][A-Z]+)[ED]$/ || $peptide =~ /^[A-Z][ED][A-Z]+\*$/){
	#if($peptide =~ /[A-Z]+QTGG[A-Z]$/){
	    print substr($peptide,1, length($peptide)-2)."\n";
	    #print "peptide is $peptide\n";
	    #print substr($peptide, 1, length($peptide)-2)."\n";
	}
    }
}
