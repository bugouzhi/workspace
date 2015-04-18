#!/usr/bin/perl
#convert raw spectrum file to mgf file
use strict;
open(RAW, $ARGV[0]) || die "cannot open raw spectrum file";
my $pmass = 0;
my $peptidePair = "";
my $charge = 0;
my %tempSpectrum = ();
my $counter =0;
while(<RAW>){
	my $line = $_;
	chomp $_;
	my @tokens = split(/\s/, $line); #we do not concatenate space as separator here because it seems in raw data some field is empty
	#$pmass = $tokens[12];
	#$peptidePair = $tokens[3];
	#$peptidePair =~ s/-/ & /;	
	if($line =~ /spectrum/){
		if((keys %tempSpectrum) > 1){
			print "BEGIN IONS\n";
			print "TITLE=spectrum_$counter.raw\n";
			print "CHARGE=$charge\n";
			#print "PEPMASS=$pmass\n";
			print "PEPMASS=".(($pmass+$charge*1.007276)/$charge)."\n";
			print "PEPSEQ=$peptidePair\n";
			foreach my $mass(sort {$a <=> $b} (keys %tempSpectrum)){
				print $mass."\t".$tempSpectrum{$mass}."\n";
			}
			print "END IONS\n\n";
		}
		$pmass = $tokens[12];
		#print "line is $line\n";
		#print "stored mass is $pmass\n";
		$peptidePair = $tokens[3];
		#print "before peptides are: $peptidePair\n";
		$peptidePair =~ s/-/ & /;
		#print "after peptides are: $peptidePair\n";
		$charge = 0; #reset charge
		%tempSpectrum = ();
		$counter++;
	}else{
		#save the spectrum line, parse the line to find out charge
		my $annotation = $tokens[0];
		$annotation =~ /_plus(\d+)/;
		#print "label is $annotation\n";
		#print "charge is $1\n";
		if($charge < $1){
			$charge = $1;
		}
		#print "key is $tokens[3]; value is $tokens[5]\n";
		$tempSpectrum{$tokens[3]} = $tokens[5];
	}
}
