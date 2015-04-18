#!/usr/bin/perl
use strict;
open(RESULT, $ARGV[0]) || die "cannot open reusult file";
open(UPS, $ARGV[1]) || die "cannot open UPS file";

my $pepInd=7;
my $protInd=8;
my $scoreInd=11;
my $UPSProt = "";

while(<UPS>){
    my $line = $_;
    chomp $line;
    if($line !~ /^>/){
	$UPSProt = $UPSProt.$line;
    }
}

while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\t/, $line);
    my $pept = @tokens[$pepInd];
    $pept =~ s/[0-9\.\+]//g;
    $pept = substr($pept, 1, length($pept)-2);
    print $pept."\n";
    if($UPSProt =~ $pept){
	@tokens[$protInd] = "UPS__".$tokens[$protInd];
    }
    print @tokens[$pepInd]."\t".@tokens[$protInd]."\n";
}
