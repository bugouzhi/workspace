#!/usr/bin/perl
#split annotation from MXDB into several part
#use to perform cross-validation

use strict;

open(RESULT, $ARGV[0]) || die "cannot open result";

my $groups = $ARGV[1];
my @charge2 = ();
my @charge3 = ();
my @charge4 = ();

while(<RESULT>){
    my $line = $_;
    my @tokens = split(/\s+/, $line);
    if($tokens[6] < 550){
	push(@charge4, $line);
    }
    
    if($tokens[6] > 550 && $tokens[6] < 800){
	push(@charge3, $line);
    
    }
    
    if($tokens[6] > 800){
	push(@charge2, $line);
    }
}

my $j = 0;
my $k = 0;
my $l = 0;
for(my $i = 1; $i <= $groups; $i++){
    open(OUT, ">annotation_part$i.txt");
    while($j < $i*$#charge2/$groups){
	print OUT $charge2[$j];
	$j++;
    }
    while($k < $i*$#charge3/$groups){
	print OUT $charge3[$k];
	$k++;
    }
    while($l < $i*$#charge4/$groups){
	print OUT $charge4[$l];
	$l++;
    }
    close(OUT);
}






