#!/usr/bin/perl
#extract a subset of results

use strict;
open(RESULTS, $ARGV[0]);
open(IDS, $ARGV[1]);

my @ids = <IDS>;
my %table = ();
my $mode = $ARGV[2];
my $index = $ARGV[3];
@table{@ids} = 1;


while(<RESULTS>){
    my $line = $_;
    my @tokens = split(/\s+/, $line);
    if(exists $table{$tokens[$index-1]."\n"} && $mode){
	print $line;
    }elsif(!exists $table{$tokens[$index-1]."\n"} && $mode == 0){
	print $line;
    }
    
}


