#!/usr/bin/perl
use strict;
#compare results from msplit-db 
#to see if different strategy for searching is consistent
#particularlly if we introduce filters in the pipeline, what is
#the frequencey that the best result matched to full db searches
open(RESULTREF, $ARGV[0]) || die "cannot open file\n";
open(RESULTTAR, $ARGV[1]) || die "cannot open file\n";
my %refTable = ();

while(<RESULTREF>){
	my $line = $_;
	chomp $line;
	if($line =~ /best/){
		my @tokens = split(/\s+/, $line);
		#print "key is $tokens[3] value is $tokens[7]\n";
		$refTable{$tokens[3]} = $tokens[7];
	}
}

my $count=0;
my $matchCount=0;
while(<RESULTTAR>){
	my $line = $_;
	chomp $line;
	if($line =~ /best/){
		my @tokens = split(/\s+/, $line);
		#print "lookup key is $tokens[3]\n";
		if(exists $refTable{$tokens[3]} && $tokens[8] > 30 && $tokens[10] > 0.5){
			$count++;
			print "spectrum $tokens[3] matched to $tokens[7] and $refTable{$tokens[3]}\n";
			if($tokens[7] =~ $refTable{$tokens[3]}){
				$matchCount++;
			}
		}
	}
}
print "matched frequency: ".($matchCount/$count)."\n";