#!/usr/bin/perl
#generate training data for svm-light from search output data
use strict;

open(RESULT, $ARGV[0]) || die "cannot open result file\n";
my $beginIndex=14;
my $ans1 = 2;
my $ans2 = 4;
my $annot1 = 9;
my $annot2 = 11;
while(<RESULT>){
    my $line =$_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    my $label=0;
    if($tokens[$ans1] =~ $tokens[$annot1] && $tokens[$ans2] =~ $tokens[$annot2] ||
       $tokens[$ans1] =~ $tokens[$annot2] && $tokens[$ans2] =~ $tokens[$annot1]){
	$label=1;
    }else{
	$label=2;
    }
    print $label." ";
    print "qid:$tokens[1] ";
    for(my $i = $beginIndex; $i <= $#tokens-2; $i++){
	print "".($i-$beginIndex+1).":".$tokens[$i]." ";
    }
    print "\n";
}
