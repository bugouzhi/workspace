#!/usr/bin/perl
#extract spectraST results
use strict;

open(RESULT, $ARGV[0]) || die "cannot open spectraST resutl file";

while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\t/, $line);
    my @tokens2 = split(/\./, $tokens[0]);
    $tokens[2] =~ s/\//\./g;
    if($tokens[16] !~ /DECOY/){
	print "Spectrum: Scan Number: ".($tokens2[1]+0)."\thas best match\t".$tokens[2]."\t".$tokens[11]."\n";
    }
}
