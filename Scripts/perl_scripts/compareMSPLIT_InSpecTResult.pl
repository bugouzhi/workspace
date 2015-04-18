#!/usr/bin/perl
use strict;
#compare results from msplit-db 
#to see if different strategy for searching is consistent
#particularlly if we introduce filters in the pipeline, what is
#the frequencey that the best result matched to full db searches
open(INSPECT, $ARGV[0]) || die "cannot open file\n";
open(MSPLIT, $ARGV[1]) || die "cannot open file\n";
my %refTable = ();

while(<INSPECT>){
	my $line = $_;
	chomp $line;
	my @tokens = split(/\t/, $line);
	#print "key is $tokens[1] value is $tokens[2]\n";
	my $pep2 = $tokens[2];
	$pep2 =~ /\.(\w+)\./;
	$pep2 = $1;
	$refTable{$tokens[1]} = $pep2."\.".$tokens[4];

}

my $count=0;
my $matchCount=0;
while(<MSPLIT>){
	my $line = $_;
	chomp $line;
	if($line =~ /best/){
		my @tokens = split(/\s+/, $line);
		#print "lookup key is $tokens[3]\n";
		if(exists $refTable{$tokens[3]}){
			$count++;
			my $pep = $tokens[7];
			#$pep =~ /(\w+)\.\d/;
			#$pep = $1;
			my $pep2 = $refTable{$tokens[3]};
			#$pep2 =~ /\.(\w+)\./;
			#$pep2 = $1;
			my $pep3 = $tokens[9];
			#$pep3 =~ /(\w+)\.\d/;
			#$pep3 = $1;
			print "spectrum $tokens[3] matched to $pep & $pep3 and $pep2\n";
			if($pep =~ $pep2 || $pep =~ $pep2){
				$matchCount++;
			}
		}
	}
}
print "matched frequency: ".($matchCount/$count)."($matchCount/$count)\n";
