#!/usr/bin/perl
#generate protein sequence variant according to seq variation info
#which is tab delimited fiel of position and varaintions

open(PROTEIN, $ARGV[0]) || die "cannot open protein sequence file\n";
open(SEQVAR, $ARGV[1]) || die "cannot open seq variation file\n";



