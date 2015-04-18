#!/usr/bin/perl
#get protein spectral count
use strict;
open(RESULT, $ARGV[0]) || die "cannot open result";
my $BAIT = $ARGV[1];
my $prevLine = "@";
my $rep=0;
my $prevProt="";
my $count=0;
while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\t/, $line);
    #print $line;
    if($tokens[6] eq $prevProt){
	$count+= $tokens[8];
    }else{
	my $ind = index($prevProt, "RepID=");
	if($ind <= 0){
	    $ind = index($prevProt, "|");
	}
	$prevProt = substr($prevProt, $ind+6);
	if($prevProt ne ""){
	    print $BAIT."\t".$BAIT."rep".$rep."\t".$prevProt."\t".$count."\n";
	}
	$count=$tokens[8];
    }    
    $prevProt = $tokens[6];
    #print "previous $prevProt\n";
    if($tokens[1] ne $prevLine){
	$rep++;
    }
    #print "previous line $prevLine\n";
    $prevLine = $tokens[1];
    
    
}
