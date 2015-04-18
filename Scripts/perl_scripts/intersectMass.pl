#!/usr/bin/perl
open(LIST1, $ARGV[0]);
open(LIST2, $ARGV[1]);

my @items1 = <LIST1>;
my @items2 = <LIST2>;
my %itemsTable =();
@itemsTable{@items1}=1;
foreach my $token(@items2){
    if(exists $itemsTable{$token}){
	print "found $token";
    }
    #if(!exists $itemsTable{$token}){
    #	print "not-found $token";
    #}
}


