#!/usr/bin/perl
#filter results from output file of search_and_bound method in SpectrumLib class
$resultFile = $ARGV[0];
@resultlines = `more $resultFile | grep 'best an'`;
foreach my $line (@resultlines){
	chomp $line;
	#print $line."\n";
	my $classification=0;
	my @tokens = split(/\s+/, $line);
	my $passed1 = ($tokens[12] >= 0.54 || $tokens[13] >= 0.54);
	my $diff1 = $tokens[11] - $tokens[12];	
	my $diff2 = $tokens[11] - $tokens[13];
	my $diff = ($diff1 < $diff2 ? $diff1 : $diff2);
	my $passed2 =  $passed1 && $tokens[14] <= 0.5; #make sure two candidate answer is not a duplicated
	$passed2 = $passed1 && $diff > ((($tokens[16]+0.000001))*0.0992 + 0.0115)*1.2; #make sure delta cosine is greater than a certain threshold
	$passed2 = $passed2 && $diff > ((($tokens[15]+0.000001))*0.0992 + 0.0115)*1.2; #make sure delta cosine is greater than a certain threshold
	#$passed2 = $passed2 && $tokens[11] > 0.72;
	$passed2 = $passed2 && $tokens[18] < 0.45;
	$passed2 = $passed2 && $tokens[19] < 0.45;
	$passed2 = $passed2 && $tokens[20] < 0.45;
	$passed = $passed && $diff > 0.02; #fixed threshold
	$thr = (($tokens[15]+0.000001))*0.0992 + 0.0115;
	$line = substr($line, 0, length($line)-1);
	if($passed2){
		$classification=2;
	}elsif($passed1){
		$classification=1;
	}
	print $line."\t".$thr."\t".$classification."\n";
}
