#!/usr/bin/perl
use strict;

my $dataDir = $ARGV[0];
my $outFile = $ARGV[1];

my @tokens = `ls $dataDir/*.dta`;
open(OUT, ">$ARGV[1]") || die "cannot open output file\n";

my $counter=1;
foreach my $file(@tokens){
    #print $file;
    open(DTA, $file) || die "cannot open file $file\n";
    print OUT "BEGIN IONS";
    if($file =~ /\/([^\/]+$)/){
	my $name = $1;
	#print "name is $name";
	print OUT "TITLE=$name";
	my $firstLine = <DTA>;
	my @tokens = split(/\s+/, $firstLine);
	my $charge = $tokens[1];
	my $parentmass = ($tokens[0] + ($charge-1)*1.007276)/$charge;
	print OUT "PEPMASS=$parentmass\n";
	print OUT "CHARGE=$charge\n";
	print OUT "SCAN=$counter\n";
	while(<DTA>){
	    print OUT $_;
	}
	print OUT "END IONS\n";
	$counter++;
    }
}
