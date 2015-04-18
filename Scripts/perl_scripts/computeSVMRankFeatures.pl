#!/usr/bin/perl
#re-compute svm features so we can set it up as a ranking problem
#for svm to learn
use strict;

open(RESULT, $ARGV[0]) || die "cannot open result file\n";

my @results = ();
my $currentScan = 0;
while(<RESULT>){
    my $line = $_;
    chomp $line;
    #$line = substr($line, 1, length($line)-1);
    my @tokens = split(/\s+/, $line);
    if($tokens[1] == $currentScan){
	push(@results, \@tokens);
    }else{
	my @average = ();
	if($#results > 0){
	    my @currToken1 = @{$results[0]};
	    #print "result has size $#results\n";
	    for(my $i = 15; $i <= $#currToken1; $i++){
		push(@average, $currToken1[$i]);	
	    }
	    #print "\n";
	    #print "finish first line\n";
	    
	    for(my $i = 1; $i <= $#results; $i++){
		my @currTokens = @{$results[$i]};
		for(my $j = 15; $j <=$#currTokens; $j++){
		    $average[$j-15] += $currTokens[$j];
		    #print "$currTokens[$j]\t";
		}
		#print "\n";
	    }
	    ##print "\n";
	    
	    for(my $i = 0; $i <= $#average; $i++){
		#print "$average[$i]\t";
		$average[$i] = $average[$i]/($#results+1);
		#print $average[$i]."\t";
	    }
	    #print "\n";
	    for(my $i = 0; $i <= $#results; $i++){
		my @currTokens = @{$results[$i]};
		for(my $j = 0; $j < 15; $j++){
		    print $currTokens[$j]."\t";
		}
		for(my $j = 15; $j <=$#currTokens; $j++){
		    print $currTokens[$j]."\t";
		    #print "$currTokens[$j]\t";
		}
		for(my $j = 15; $j <=$#currTokens; $j++){
		    if($average[$j-15] == 0){
			$average[$j-15] = 0.0000001;
		    }
		    print $currTokens[$j]/$average[$j-15]."\t";
		    #print "$currTokens[$j]\t";
		}
		print "\n";
	    }
        }
	@results=();
	push(@results, \@tokens);
	$currentScan = $tokens[1];
	#print "current scan $currentScan\n";
      
    }
}
