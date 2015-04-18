#!/usr/bin/perl
open(LIST1, $ARGV[0]);
open(LIST2, $ARGV[1]);

my @items1 = <LIST1>;
my @items2 = <LIST2>;
my $Index = $ARGV[2];
my %itemsTable =();
@itemsTable{@items1}=1;

my $sep = "\t";
foreach my $line(@items1){
    chomp $line;
    my @tokens = split($sep, $line);
    #print $tokens[$Index]."\n";
    chomp $tokens[$Index];
    @itemsTable{$tokens[$Index]} = 1;
}

foreach my $line(@items2){
    my @tokens = split($sep, $line);
    if(exists $itemsTable{$tokens[$Index]}){
	if($tokens[$Index] =~ /\n/){
	    print "".$tokens[$Index];
	}else{
	    print "".$tokens[$Index]."\n";
	}
    }
    #if(!exists $itemsTable{$token}){
    #	print "not-found $token";
    #}
}


