#!/usr/bin/perl
#post-process output from MSPLIT and mapeed out decoy peptides
use strict;
open(HITS, $ARGV[0]) || die "cannot open forward db file";
open(DECOY, $ARGV[1]) || die "cannot open reverse db file";
open(MSPLIT, $ARGV[2]) || die "cannot open MSPLIT result file";

my %forwardHitsTable = ();
my %reverseHitsTable = ();
my @hits = <HITS>;
@forwardHitsTable{@hits} = 1;
my @decoy = <DECOY>;
@reverseHitsTable{@decoy} = 1;
my $count = 0;
while(<MSPLIT>){
	my $line = $_;
	chomp $line;
	if($line =~ /best/){
#		if($count < 2){
#			$count++;
#			next;
#		}
		$count = 0;
		#print $line."\n";
		my @tokens = split(/\s+/, $line);
		#my $peptide1 = substr($tokens[7], 0, length($tokens[7])-2)."\n";
		#my $peptide2 = substr($tokens[9], 0, length($tokens[9])-2)."\n";
		my $peptide1 = $tokens[10]."\n";
		my $peptide2 = $tokens[12]."\n";
		$peptide1 =~ s/[0-9\.\+]//g;
		$peptide2 =~ s/[0-9\.\+]//g;
		#print "peptide one is $peptide1".(length($peptide1))."\n";
		#print "peptide two is $peptide2".(length($peptide2))."\n";
		#print "peptide one is $tokens[10] ".(length($peptide1))."\n";
		#print "peptide two is $tokens[12] ".(length($peptide2))."\n";
		
		if(!exists $forwardHitsTable{$peptide1} && exists $reverseHitsTable{$peptide1}){
			$line =~ s/$tokens[10]/r$tokens[10]/;
			#print "found $line\n";
		}	
		if(!exists $forwardHitsTable{$peptide2} && exists $reverseHitsTable{$peptide2}){
			$line =~ s/$tokens[12]/r$tokens[12]/;
			#print "found $line\n";
		}				
		print $line."\n";
	}
}
