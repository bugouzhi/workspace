#!/usr/bin/perl
use strict;
#compare results from msplit-db 
#to see if different strategy for searching is consistent
#particularlly if we introduce filters in the pipeline, what is
#the frequencey that the best result matched to full db searches
open(INSPECT, $ARGV[0]) || die "cannot open file\n";
open(MSPLITDB, $ARGV[1]) || die "cannot open file\n";
my %refTable = ();

while(<INSPECT>){
	my $line = $_;
	chomp $line;
	my @tokens = split(/\s+/, $line);
	#print "key is $tokens[1] value is $tokens[2]\n";
	$refTable{$tokens[1]} = $tokens[2];

}

my $count=0;
my $matchCount=0;
while(<MSPLITDB>){
	my $line = $_;
	chomp $line;
	if($line =~ /best/){
		my @tokens = split(/\s+/, $line);
		#print "lookup key is $tokens[3]\n";
		if(exists $refTable{$tokens[3]}){
			$count++;
			my $pep = $tokens[7];
			$pep =~ /(\w+)\.\d/;
			$pep = $1;
			my $pep1 = $tokens[9];
			$pep1 =~ /(\w+)\.\d/;
			$pep1 =$1;
			my $pep2 = $refTable{$tokens[3]};
			$pep2 =~ /\.(\w+)\./;
			$pep2 = $1;
			print "spectrum $tokens[3] matched to $pep and $pep2\n";
			if($pep =~ $pep2 || $pep1 =~ $pep2){
				$matchCount++;
			}
		}
	}
}
print "matched frequency: ".($matchCount/$count)."($matchCount/$count)\n";
