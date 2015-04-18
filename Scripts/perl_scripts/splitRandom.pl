#!/usr/bin/perl
use strict;
#split mixdb results randomly into parts
#use to perform cross-validation

open(RESULT, $ARGV[0]) || die "cannot open result file\n";
my $partCount = $ARGV[1];

my @fhandles = ();
for(my $i = 1; $i <= $partCount; $i++){
    open(my $fh, ">out_part$i.txt") || die "cannot open output file\n";
    push(@fhandles, $fh);
}

my $prevScan = -1;
my $bin = -1;
while(<RESULT>){
    my $line = $_;
    my @tokens = split(/\s+/, $line);
    if($tokens[1] == $prevScan){
	my $fh = $fhandles[$bin];
	print $fh $line;
    }else{
	$bin = int(rand($partCount));
	#print "bin: $bin\n";
	my $fh = $fhandles[$bin];
	print $fh $line;
	$prevScan = $tokens[1];
    }
}
