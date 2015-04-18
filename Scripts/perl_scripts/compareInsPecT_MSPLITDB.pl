#!/usr/bin/perl
use strict;
#compare results from msplit-db 
#to see if different strategy for searching is consistent
#particularlly if we introduce filters in the pipeline, what is
#the frequencey that the best result matched to full db searches
open(INSPECT, $ARGV[0]) || die "cannot open file\n";
open(MSPLITDB, $ARGV[1]) || die "cannot open file\n";
my %refTable1 = ();
my %refTable2 = ();


while(<MSPLITDB>){
    my $line=$_;
    chomp $line;
   if($line =~ /best/){
       my @tokens = split(/\s+/, $line);
       my $pep = $tokens[7];
       $pep =~ /(\w+)\.\d/;
       $pep = $1;
       $refTable1{$tokens[3]} = $pep;
       my $pep1 = $tokens[9];
       $pep1 =~ /(\w+)\.\d/;
       $pep1 =$1;
       $refTable2{$tokens[3]} = $pep1;
   } 
}

my $count=0;
my $matchCount=0;

while(<INSPECT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    #print "key is $tokens[1] value is $tokens[2]\n";
    if(exists $refTable1{$tokens[1]}){
	$count++;
	my $pep1 = $refTable1{$tokens[1]};
	my $pep2 = $refTable2{$tokens[1]};
	my $pep = $tokens[2];
	$pep =~ /\.(\w+)\./;
	$pep = $1;
	print "spectrum $tokens[1] matched to $pep and $pep1\n";
	if($pep =~ $pep1 || $pep =~ $pep1){
	    $matchCount++;
	}
    }
}

print "matched frequency: ".($matchCount/$count)."($matchCount/$count)\n";
