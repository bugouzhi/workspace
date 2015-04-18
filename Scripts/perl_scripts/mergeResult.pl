#!/usr/bin/perl
#merge result from annotation and parentmass check
use strict;

open(RESULT1, $ARGV[0]);
open(RESULT2, $ARGV[1]);

my %table = ();

while(<RESULT1>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    $table{$tokens[0]}=$line;
}

while(<RESULT2>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    if(exists $table{$tokens[0]}){	
		print substr($line, 0,length($line)-2)."\t".$table{$tokens[0]}."\n";
    }else{
		print substr($line, 0,length($line)-2)."\t"."not found\n";
    }
}
