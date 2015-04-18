#!/usr/bin/perl
use strict;
#count number of misscleavage

open(RESULT, $ARGV[0]) || die "cannot open result file\n";
my $pepIndex = $ARGV[1];

while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\t/, $line);
    my $peptide = $tokens[$pepIndex];
    if($peptide =~ /\w\.(.+)\.\w/){
	$peptide = $1;
	my $count=0;
	while($peptide =~ /[K]/g){
	    $count++;
	}
	#while($peptide =~ /[RK]/g){
	#    $count++;
	#}
	print $line."\t".$count."\n";
    }
}
