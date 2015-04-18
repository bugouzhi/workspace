#!/usr/bin/perl
use strict;
open(FILE, $ARGV[0]) || die "cannot open output file";
my %table=();

while(<FILE>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/[\s_]/, $line);
    #print "key is $tokens[1]\n";
    if($tokens[5] < 0.0103){
	if(exists $table{$tokens[1]}){
	    my @items = @{$table{$tokens[1]}};
	    push(@items, $line);
	    $table{$tokens[1]} = \@items;
	}else{
	    my @items = ();
	    push(@items, $line);
	    $table{$tokens[1]} = \@items;
	}
    }

}

#print "table has entries ".%table."\n";

foreach my $key (keys %table){
    my @items = @{$table{$key}};
    if($#items > 0){
	print $key;
	foreach my $item(@items){
	    #print "\t".substr($item, 0, length($item)-1);
	    my @tokens = split(/\s+/, $item);
	    my @tokens2 = split(/_/, $tokens[0]);
	    @tokens = split(/\./, $tokens[3]);
	    print "\t".$tokens[1].".".$tokens2[2];
	}
	print "\n";
    }
}
