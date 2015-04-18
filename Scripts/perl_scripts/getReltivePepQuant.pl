#!/usr/bin/perl
#compute relative peptide abundance for proteins in UPS1 vs UPS2

use strict;


open(TABLE, $ARGV[0]) || die "cannot open file";

my $prot="A";
my @abundance1=();
my @abundance2=();
my $abundInd1=1;
my $abundInd2=2;
my $protInd=5;
<TABLE>;
while(<TABLE>){
    my @tokens = split(/\t/, $_);
    #print $tokens[5]."\n";
    #print $prot."\n";
    my $sameProb = 0.8;
    if(rand() < $sameProb){
	#$tokens[$protInd] = 0;
    }else{
	#$tokens[$protInd] = 1;
    }
    if($tokens[$protInd] eq $prot){
	push(@abundance1, $tokens[$abundInd1]);
	push(@abundance2, $tokens[$abundInd2]);
    }else{
	#print "here\n";
	#print "$#abundance1\n";
	for(my $i = 0; $i <= $#abundance1;  $i++){
	    for(my $j = $i+1; $j <= $#abundance1; $j++){
		my $ratio1 = $abundance1[$i] / $abundance1[$j];
		my $ratio2 = $abundance2[$i] / $abundance2[$j];
		print $tokens[0]."\trelative ratio:\t".$ratio1."\t".$ratio2."\t".$abundance1[$i]."\t".$abundance1[$j]."\n";
	    } 
	}
	$prot = $tokens[$protInd];
	#print $prot."\n";
	@abundance1=();
	@abundance2=();
	push(@abundance1, $tokens[$abundInd1]);
	push(@abundance2, $tokens[$abundInd2]);
    }
}
