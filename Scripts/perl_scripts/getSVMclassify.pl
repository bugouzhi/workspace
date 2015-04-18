#!/usr/bin/perl
#generate training data for svm-light from search output data
use strict;

my $svminput="temp_svmin.txt";
my $svmresult="temp_svmout.txt";

open(RESULT, $ARGV[0]) || die "cannot open result file\n";
open(SVMIN, $svminput) || die "cannot open svm temp input file\n";

my $beginIndex=15;
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
    for(my $i = $beginIndex; $i <= $#tokens-2; $i++){
	print SVMIN "".($i-$beginIndex+1).":".$tokens[$i]." ";
    }
    print "\n";
}
close(RESULT);
close(SVMIN);

`./svm-classify.exe temp_svmin.txt $svmmodel $svmout`;

open(RESULT, $ARGV[0]) || die "cannot open result file-2\n";
open(SVMRESULT, $svmresult) || die "cannot open svm-result file\n";
while(<RESULT>){
    my $line = $_;
    my $svm = <SVMRESULT>;
    chomp $line;
    chomp $svm;
    print $line."\t"$svm."\n";
    
}
