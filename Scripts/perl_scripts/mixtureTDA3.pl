#!/usr/bin/perl
#we perform TDA analysis for peptide-peptide-spectrum match: ppsm
#note in this version we only consider correct ppsms as those where both are correct
use strict;

open(RESULT, $ARGV[0]) || die "cannot open results";
my $FDR = $ARGV[1];
my $mode= $ARGV[2];
my @target1=();
my @target2=();
my @decoy1=();
my @decoy2=();
my @decoy22=();
my $scoreIndex1=33;
my $scoreIndex2=33;
my $threshold1=10000;
my $threshold2=10000;
my $pepIndex1 = 10;
my $pepIndex2 = 12;


while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    if($tokens[$pepIndex1] =~ /[rxX]/ && $tokens[$pepIndex2] =~ /[rxX]/){
	#print $tokens[$scoreIndex1]."\n";
	push(@decoy1, $mode*$tokens[$scoreIndex1]);
    }else{
	push(@target1, $mode*$tokens[$scoreIndex1]);
    }
}

close(RESULT);
$threshold1 = getThreshold(\@target1, \@decoy1, 1.0);#$FDR/3);  #no first stage at this moment
print "threshold1: $threshold1\n";  

open(RESULT, $ARGV[0]) || die "cannot open results";
while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    if($mode*$tokens[$scoreIndex1] > $threshold1 && $tokens[15] > -100){
	if($tokens[$pepIndex1] =~ /[rxX]/ && $tokens[$pepIndex2] !~ /[rxX]/ || $tokens[$pepIndex1] !~ /[rxX]/ && $tokens[$pepIndex2] =~ /[rxX]/){
	    print tokens[$scoreIndex1]."\n";
	    push(@decoy2, $mode*$tokens[$scoreIndex2]);
	}elsif($tokens[$pepIndex1] !~ /[rxX]/ && $tokens[$pepIndex2] !~ /[rxX]/){
	    push(@target2, $mode*$tokens[$scoreIndex2]);
	}else{
	    push(@decoy22, $mode*$tokens[$scoreIndex2]);
	}
    }
}
close(RESULT);
$threshold2 = getThreshold(\@target2, \@decoy2, $FDR);
$threshold2 = getThreshold12(\@target2, \@decoy2, \@decoy22, $FDR);


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


sub getThreshold12{
    
    my @sortedTarget = sort {$a <=> $b} @{$_[0]}; 
    my @sortedDecoy = sort {$a <=> $b} @{$_[1]};
    my @sortedDecoy2 = sort {$a <=> $b} @{$_[2]};
    my $fdr = $_[3];
    print "number of targets: $#sortedTarget\n";
    print "number of decoys: ".($#sortedDecoy+$#sortedDecoy2)."\n";
    my $i=0;
    my $j=0;
    my $k=0;
    my $threshold=10000;
    my $targetCount = 0;
    my $decoyCount=0;
    my $decoyCount2=0;
    while($i <= $#sortedTarget && $j <= $#sortedDecoy){
	#print "$i\t$j\t$k\n";
	#print $sortedTarget[$i]."\t".$sortedDecoy[$j]."\n";
        if($sortedTarget[$i] > $sortedDecoy[$j]){
	    $j++;
	    next;
	}
	
	while($sortedDecoy2[$k] <= $sortedDecoy[$j+1] && $k <= $#sortedDecoy2){
	    $k++;
	}
		
	#print "ratio: ".(($#sortedDecoy + 1 - $j) / ($#sortedTarget + 1 -$i))."\n";
        $decoyCount =$#sortedDecoy + 1 - $j;
	$decoyCount2 = $#sortedDecoy2 + 1 - $k;
	$targetCount = $#sortedTarget + 1 -$i;
	
	#print "count: $targetCount\t$decoyCount\t$decoyCount2\n"; 
	if((($decoyCount - $decoyCount2) / $targetCount) < $fdr && ($decoyCount2 / $targetCount) < ($fdr*3)){ #make sure error is not far off
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
    print "counts: T=$targetCount\tD=$decoyCount\tD2=$decoyCount2\n"; 
    print "number of accepted ids: ".($#sortedTarget+1 - $i)."\tthreshold is $threshold\n";
    return $threshold;
}
