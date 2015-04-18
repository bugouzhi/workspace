#!/usr/bin/perl
open(MGF, $ARGV[0]);
open(OUT, ">$ARGV[1]");
while(<MGF>){
	my $line = $_;
	chomp $line;
	if($line =~ /\d/){
		my @tokens = split(/\s+/, $line);
		print OUT "1\t$tokens[0]\t\?\t$tokens[1]\n";
	}
}