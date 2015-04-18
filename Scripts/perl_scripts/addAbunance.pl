#!/usr/bin/perl
use strict;

open(RESULT, $ARGV[0]) || die "cannot open file";
open(ABUN, $ARGV[1]) || die "cannot open file2";

my %map1 = {};
my $pepInd=4;
my $proInd=5;

while(<ABUN>){
    my $line = $_;
    chomp $line;
    if($line =~ /^>/){
	my @tokens = split(/\s+/, $line);
	my @tokens2 = split(/[|_]/, $tokens[0]);
	#print "key-map: ".substr($tokens[0],1)."\t".$tokens[$#tokens]."\n";
	#print "key-map: ".$tokens2[2]."\t".$tokens[$#tokens]."\n";
	#$map1{substr($tokens[0],1)} = $tokens[$#tokens];
	$map1{$tokens2[2]} = $tokens[$#tokens];
    }
}

while(<RESULT>){
    my $line = $_;
    chomp $line;
    $line = substr($line, 0, length($line)-2);
    my @tokens = split(/\t/, $line);
    my @tokens2 = split(/[|_]/, $tokens[$proInd]);
    my $abundance = $map1{$tokens2[2]};
    #print (exists $map1{"CAH1"})."\n";
    print "look: ".$tokens2[2]."\n";
    print $line."\t";
    print $abundance."\n";
}
