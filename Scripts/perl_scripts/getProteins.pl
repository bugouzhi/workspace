#!/usr/bin/perl
use strict;

open(IDS, $ARGV[0]) || die "cannot open protein IDs file\n";
open(PROTEINS, $ARGV[1]) || die "cannot open protein fasta file\n";

my @IDS=();
while(<IDS>){
    chomp $_;
    push(@IDS, $_);
}

my $match=0;
while(<PROTEINS>){
    my $line=$_;
    if($line =~ /^>/){
	$match=0;
	foreach my $id(@IDS){
	    #print "id is $id\n";
	    my @tokens = split(/\s+|\|/, $id);
	    $tokens[1]=$id;
	    my $name = $tokens[1];
	    $name =~ s/\|/\\\|/g;
	    #print "id: $name line ".substr($line,1)."\n";
	    if($line =~ $name){
		#print "match: id: $tokens[1] line: $line\n";
		$match=1;
	    }
	}
    }
    if($match==1){
	print $line;
    }
}
