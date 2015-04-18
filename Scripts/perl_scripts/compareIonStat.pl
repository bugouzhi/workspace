#!/usr/bin/perl
#compare ions statistics between TheoreticalSpectrum matching in m-split and inspect
open(INSPECT, $ARGV[0]);
open(MSPLIT, $ARGV[1]);

my %table = ();
while(<INSPECT>){
	my $line = $_;
	chomp $line;
	my @tokens = split(/\t/, $line);
	$table{$tokens[1]} = \@tokens;
}
close(INSPECT);

while(<MSPLIT>){
	my $line = $_;
	chomp $line;
	my @tokens = split(/\s|\(/, $line);
	#print "line: ".$line."\n";
	#print $tokens[2]."\t".$tokens[9]."\n";
	my @tokens2 = @{$table{$tokens[2]}};
	#print $tokens[2]."\t".$tokens2[11]."\n";
	print $tokens[0]."\t".$tokens[1].$tokens[2]." ions fraction difference: ".($tokens[9]-$tokens2[11])."\n"
}
