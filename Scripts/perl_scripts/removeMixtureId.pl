#!/usr/bin/perl
#remove possible mixture ids from probidtree result
use strict;
open(RESULT, $ARGV[0]) || "die cannot open result file";
my %table=();
while(<RESULT>){
    my $line = $_;
    chomp $line;
    if($line =~ /closest match/){
	my @tokens = split(/\s+/, $line);
	my $scanNum = $tokens[3];
	if(exists $table{$scanNum}){
	    my $count= $table{$scanNum};
	    $count++;
	    $table{$scanNum} = $count;
	}else{
	    $table{$scanNum} = 1;
	}
    }
}
close(RESULT);

open(RESULT, $ARGV[0]);

while(<RESULT>){
    my $line = $_;
    chomp $line;
    if($line =~ /closest match/){
	my @tokens  = split(/\s+/, $line);
	if($table{$tokens[3]} == 1){
	    #print "scan $tokens[3] appear more than once\n";
	    print $line."\n";
	}
    }
}

