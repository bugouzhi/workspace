#!/usr/bin/perl
use strict;
#compare results from msplit-db 
#to see if different strategy for searching is consistent
#particularlly if we introduce filters in the pipeline, what is
#the frequencey that the best result matched to full db searches
open(SPECTRAST, $ARGV[0]) || die "cannot open file\n";
open(MSPLIT, $ARGV[1]) || die "cannot open file\n";
my %refTable = ();

while(<SPECTRAST>){
	my $line = $_;
	chomp $line;
	my @tokens = split(/\s+/, $line);
	my @tokens2 = split(/\./, $tokens[0]);
	my @tokens3 = split(/\//, $tokens[2]);
	#print "key is ".($tokens2[1]+0)." value is $tokens3[0]\n";
	$refTable{($tokens2[1]+0)} = $tokens3[0];

}

my $count=0;
my $matchCount=0;
while(<MSPLIT>){
	my $line = $_;
	chomp $line;
	if($line =~ /best/){
		my @tokens = split(/\s+/, $line);
		#print "lookup key is $tokens[4]\n";
		if(exists $refTable{$tokens[4]}){
			$count++;
			my $pep = $tokens[5];
			$pep =~ /(\w+)\.\d/;
			$pep = $1;
			my $pep2 = $refTable{$tokens[4]};
			$pep2 =~ /\.(\w+)\./;
			$pep2 = $1;
			my $pep3 = $tokens[7];
			$pep3 =~ /(\w+)\.\d/;
			$pep3 = $1;
			print "spectrum $tokens[4] matched to $pep & $pep3 and $pep2\n";
			if($pep =~ $pep2 || $pep3 =~ $pep2){
				$matchCount++;
			}
		}
	}
}
print "matched frequency: ".($matchCount/$count)."($matchCount/$count)\n";
