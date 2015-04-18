#!/usr/bin/perl
#labelled the result from MSPLIT into no-match, single match and mixture match respectively
open(MSPLIT, $ARGV[0]) || die "cannot open MSPLIT result file";
while(<MSPLIT>){
	my $line = $_;
	chomp $line;
	my @tokens= split(/\s+/, $line);
	$line = substr($line, 0, length($line)-2);
	if($#tokens > 26){
		if($tokens[1] =~ /Scan/){
			print $line."\t1\n";
		}
		if(($tokens[1] =~ $tokens[7] && $tokens[3] =~ $tokens[9])
		|| ($tokens[1] =~ $tokens[9] && $tokens[3] =~ $tokens[7])){
			print $line."\t3\n";
		}elsif($tokens[1] =~ $tokens[7] || $tokens[3] =~ $tokens[7]){
			print $line."\t2\n";
		}else{
			print $line."\t1\n";
		}
	}else{
		if($tokens[1] =~ $tokens[5] || $tokens[1] =~ $tokens[7]){
			print $line."\t2\n";
		}else{
			print $line."\t1\n";
		}
	}
}