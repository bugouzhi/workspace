#!/usr/bin/perl
#compare result from different search
#line up the results one after another by scan number

use strict;

open(RESULT1, $ARGV[1]) || "cannot open result file1";
open(RESULT2, $ARGV[0]) || "cannot open result file2";
my $direction1 = $ARGV[4];
my $direction2 = $ARGV[5];
my $threshold1 = $direction1*$ARGV[2];
my $threshold2 = $direction2*$ARGV[3];


my %table = ();
my $key=3;
my $Id=7;
my $score=8;

while(<RESULT1>){
    my $line = $_;
    #print "line is $line\n";
    my @tokens = split(/\s+/, $line);
    #print "key is $tokens[$key]\n";
    $table{$tokens[$key]} = $line;
}

my @index=(7,8);
my $overlapSame=0;
my $overlapDiff=0;
my $diffSame1=0;
my $diffDiff1=0;
my $diffSame2=0;
my $diffDiff2=0;

while(<RESULT2>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    if(exists $table{$tokens[$key]}){
	@tokens = split(/\s+/, $line);
	print $tokens[$key]."\t";
	for(my $i=0; $i <= $#index; $i++){
	    print $tokens[$index[$i]]."\t";
	}
	print "\t";
	my $line2 = $table{$tokens[$key]};
	my @tokens2 = split(/\s+/, $line2);
	for(my $i=0; $i <= $#index; $i++){
	    print $tokens2[$index[$i]]."\t";
	}
	print "\n";
	#not support PTM now
	$tokens[$Id] =~ s/[\[\]\+\-\.0-9]//g;
	$tokens2[$Id] =~ s/[\[\]\+\-\.0-9]//g;
	
	my $score1 = $direction1*$tokens[$score];
	my $score2 = $direction2*$tokens2[$score];
	#print "$score1\t$threshold1\t$score2\t$threshold2\n";
	if($score1 > $threshold1 && $score2 > $threshold2){
	    if($tokens[$Id] =~ $tokens2[$Id]){
		$overlapSame++;
	    }else{
		$overlapDiff++;
	    }
	}elsif($score1 > $threshold1){
	    if($tokens[$Id] =~ $tokens2[$Id]){
		$diffSame1++;
	    }else{
		$diffDiff1++;
	    }
	}elsif($score2 > $threshold2){
	    if($tokens[$Id] =~ $tokens2[$Id]){
		$diffSame2++;
	    }else{
		$diffDiff2++;
	    }
	}

	delete $table{$tokens[$key]};
    }else{
	print $tokens[$key]."\t";
	for(my $i=0; $i <= $#index; $i++){
	    print $tokens[$index[$i]]."\t";
	}
	print "\t";
	for(my $i=0; $i <= $#index; $i++){
	    print "missing\t";
	}
	print "\n";
    }
}

#printing remaining
foreach my $key (keys %table){
    print $key."\t";
    my $line2 = $table{$key};
    my @tokens = split(/\s+/, $line2);
    for(my $i=0; $i <= $#index; $i++){
	print "missing\t";
    }
    print "\t";
    for(my $i=0; $i <= $#index; $i++){
	print $tokens[$index[$i]]."\t";
    }
    print "\n";
}

print "Summary: $overlapSame\t$overlapDiff\t$diffSame1\t$diffDiff1\t$diffSame2\t$diffDiff2\n";
