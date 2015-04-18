#!/usr/bin/perl
use strict;
#compare results from msplit-db 
#to see if different strategy for searching is consistent
#particularlly if we introduce filters in the pipeline, what is
#the frequencey that the best result matched to full db searches
open(MSPLIT, $ARGV[0]) || die "cannot open file\n";
open(MSPLITDB, $ARGV[1]) || die "cannot open file\n";
my %refTable = ();

while(<MSPLIT>){
	my $line = $_;
	if($line =~ /best/){
	    chomp $line;
	    my @tokens = split(/\s+/, $line);
	    #print "key is $tokens[3] value is $tokens[7]\n";
	    #if($tokens[11] > -1000){
	    $refTable{$tokens[3]} = $tokens[7]." & ".$tokens[9];
	    #}
	}
}

my $count=0;
my $matchCount=0;
while(<MSPLITDB>){
	my $line = $_;
	chomp $line;
	if($line =~ /best/){
		my @tokens = split(/\s+/, $line);
		#print "lookup key is $tokens[3]\t".length($tokens[3])."\n";
		if(exists $refTable{$tokens[3]}){
			$count++;
			my $pep = $tokens[7];
			my $pep2 = $tokens[9];
			my $match = $refTable{$tokens[3]};
			my @peps = split(/ & /, $match);
			#print "spectrum $tokens[3] matched to $pep & $pep2"." \t  "."$match\n";
			if($pep =~ $peps[0] || $pep =~ $peps[0] || $pep =~ $peps[0] || $pep =~ $peps[0]){
			    #print "spectrum $tokens[3] matched to $pep & $pep2"." \t  "."$match\n";
			    $matchCount++;
			}else{
			    #print "spectrum $tokens[3] differentially matched to $pep & $pep2"." \t  "."$match\n";
			}
			if($pep =~ $peps[0] && $pep2 =~ $peps[1] || $pep =~ $peps[1] && $pep2 =~ $peps[0]){
			    print "spectrum $tokens[3] equivalently matched to $pep & $pep2"." \t  "."$match\n";
			}else{
			    print "spectrum $tokens[3] differentially matched to $pep & $pep2"." \t  "."$match\n";
			}
		    }else{
			my $pep = $tokens[7];
			my $pep2 = $tokens[9];
			#print "spectrum $tokens[3] differentially matched to $pep & $pep2\tnothing\n";
		    }
	}
}
print "matched frequency: ".($matchCount/$count)."\t$matchCount/$count\n";
