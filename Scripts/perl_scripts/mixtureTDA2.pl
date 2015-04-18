#!/usr/bin/perl
#we perform TDA analysis for mixture spectrum
use strict;

open(RESULT, $ARGV[0]) || die "cannot open results";
my $FDR = $ARGV[1];
my $mode= $ARGV[2];
my @target1=();
my @target2=();
my @decoy1=();
my @decoy2=();
my $scoreIndex1=40;
my $scoreIndex2=39;
my $threshold1=10000;
my $threshold2=10000;
my $pepIndex1 = 12;
my $pepIndex2 = 14;



while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    #$tokens[$scoreIndex1]=$tokens[$scoreIndex1]*$tokens[$scoreIndex2];
    if($tokens[$pepIndex1] =~ /^[rxX]/ && $tokens[$pepIndex2] =~ /^[rxX]/){
	#print $tokens[$scoreIndex1]."\n";
	push(@decoy1, $mode*$tokens[$scoreIndex1]);
    }elsif($tokens[$pepIndex1] !~ /^[rxX]/ && $tokens[$pepIndex2] !~ /^[rxX]/){
	push(@target1, $mode*$tokens[$scoreIndex1]);
    }
}

close(RESULT);
$threshold1 = getThreshold(\@target1, \@decoy1, $FDR);
print "threshold1: $threshold1\n";

open(RESULT, $ARGV[0]) || die "cannot open results";
while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    #$tokens[$scoreIndex1]=$tokens[$scoreIndex1]*$tokens[$scoreIndex2];
    if($mode*$tokens[$scoreIndex1] > $threshold1 && $tokens[15] > -100){
	if($tokens[$pepIndex1] =~ /^[rxX]/ && $tokens[$pepIndex2] !~ /^[rxX]/ || $tokens[$pepIndex1] !~ /^[rxX]/ && $tokens[$pepIndex2] =~ /^[rxX]/){
	    print tokens[$scoreIndex1]."\n";
	    push(@decoy2, $mode*$tokens[$scoreIndex2]);
	}elsif($tokens[$pepIndex1]!~ /^[rxX]/ && $tokens[$pepIndex2] !~ /^[rxX]/){
	    push(@target2, $mode*$tokens[$scoreIndex2]);
	}
    }
}
close(RESULT);
$threshold2 = getThreshold(\@target2, \@decoy2, $FDR*2);


sub getThreshold{
    
    my @sortedTarget = sort {$a <=> $b} @{$_[0]}; 
    my @sortedDecoy = sort {$a <=> $b} @{$_[1]};
    my $fdr = $_[2];
    print "number of targets: $#sortedTarget\n";
    print "number of decoys: $#sortedDecoy\n";
    my $i=0;
    my $j=0;                                                    
    my $threshold=10000;
  
    while($i <= $#sortedTarget && $j <= $#sortedDecoy){
	#print $i."\t".$j."\n";
	#print $sortedTarget[$i]."\t".$sortedDecoy[$j]."\n";
        if($sortedTarget[$i] > $sortedDecoy[$j]){
	    $j++;
	    next;
	}
	#print "target: ".$sortedTarget[$i]."\n";
	#print "ratio: ".(($#sortedDecoy + 1 - $j) / ($#sortedTarget + 1 -$i))."\n";
        if(($#sortedDecoy + 1 - $j) / ($#sortedTarget + 1 -$i) < $fdr){
	    $threshold = $sortedTarget[$i];
	    last;
	}         
	if($sortedTarget[$i] <= $sortedDecoy[$j]){
	    $i++;
	    next;
	}
    }
    if($j >= $#sortedDecoy){
        $threshold = $sortedTarget[$i];
    }
    print "number of accepted ids: ".($#sortedTarget+1 - $i)."\tthreshold is $threshold\n";
    return $threshold;
}


sub getThreshold2{
    
    my @sortedTarget = sort {$a <=> $b} @{$_[0]}; 
    my @sortedDecoy = sort {$a <=> $b} @{$_[1]};
    my $fdr = $_[2];
    print "number of targets: $#sortedTarget\n";
    print "number of decoys: $#sortedDecoy\n";
    my $i=0;
    my $j=0;                                                    
    my $threshold=10000;
  
    while($i <= $#sortedTarget && $j <= $#sortedDecoy){
	#print $i."\t".$j."\n";
	#print $sortedTarget[$i]."\t".$sortedDecoy[$j]."\n";
        if($sortedTarget[$i] > $sortedDecoy[$j]){
	    $j++;
	    next;
	}
	#print "ratio: ".(($#sortedDecoy + 1 - $j) / ($#sortedTarget + 1 -$i))."\n";
        if(($#sortedDecoy + 1 - $j) / ($#sortedTarget + 1 -$i) < $fdr){
	    $threshold = $sortedTarget[$i];
	    last;
	}         
	if($sortedTarget[$i] <= $sortedDecoy[$j]){
	    $i++;
	    next;
	}
    }
    if($j >= $#sortedDecoy){
        $threshold = $sortedTarget[$i];
    }
    print "number of accepted ids: ".($#sortedTarget+1 - $i)."\tthreshold is $threshold\n";
    return $threshold;
}
