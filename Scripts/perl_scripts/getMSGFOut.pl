#!/usr/bin/perl
use strict;

open(RESULT, $ARGV[0]) || die "cannot open MSGF result file\n";

while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\t/, $line);
    $tokens[3] =~ s/C\-0\.031/C/g;
    $tokens[3] =~ /\*\.(.+)\.\*/;
    my $pep = $1;
    if($tokens[5] !~ /Decoy/ && $line !~ /^\#/){
	print "Spectrum Scan Number:\t$tokens[2]\thas best match:\t$pep\.$tokens[4]\t$tokens[6]\n";
    }
}
