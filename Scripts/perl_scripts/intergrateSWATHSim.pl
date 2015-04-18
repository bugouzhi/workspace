#!/usr/bin/perl
open(RESULT, $ARGV[0]) || die "cannot open result files";

my $peptide = "";
my $protein = "";
my $total = 0;
my $mean = 0;
my $var = 0;
my @scores = ();
my $best = 0;
my $prevScan = 0;
my $beginScan = 0;
while(<RESULT>){
    my $line = $_;
    chomp $line;
    @tokens = split(/\t/, $line);
    #print "peptide is $tokens[4]\n";
    #print "scan $prevScan\n";
    if($peptide eq $tokens[4] && $tokens[1] - $prevScan < 75){
	if($tokens[7] > 0.75 && $tokens[11] > 8){
	    #print "here\n";
	    #print "$#scores\n";
	    push(@scores, $tokens[7]);
	}
    }else{
	if($#scores >= 0){
	    $total = 0;
	    foreach my $score(@scores){
		$total += $score;
		if($score > $best){
		    $best = $score;
		}
	    }
	    $mean = $total/($#scores+1);
	    foreach my $score(@scores){
		$var += ($score-$mean)*($score-$mean);
	    }
	    $var = sqrt($var)/($#scores+1);
	    print "Scan\t".$beginScan."\t".$prevScan."\tpeptide:\t".$peptide."\t".$protein."\t".$mean."\t".$total."\t".$var."\t".$best."\t".($#scores+1)."\n";	    
	}
	@scores = ();
	push(@scores, $tokens[7]);
	$best = 0;
	$peptide = $tokens[4];
	$protein = $tokens[8];
	$beginScan = $tokens[1];
      }
    $prevScan = $tokens[1];
}
