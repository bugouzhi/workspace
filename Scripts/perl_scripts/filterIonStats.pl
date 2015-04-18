#!/usr/bin/perl
open(RESULT, $ARGV[0]) || die "cannot open result files";
while(<RESULT>){
	my $line = $_;
	chomp $line;
	if($line =~ /ions fraction/){
		my @tokens = split(/\s+|\)|\(/, $line);
		my @savedLines = ();
		push(@savedLines, $line);
		#print "$line is: $line\n";
		#print "threshold is $tokens[14]\n";
		if($tokens[14] > 0.8 && $tokens[7] < 5){
			my @ionsCoverage = ();
			for(my $i = 0; $i < 4; $i++){
				$line = <RESULT>;
				chomp $line;
				push(@savedLines, $line);
				@tokens = split(/\t/, $line);
				#print "coverage is $tokens[2]\n";
				push(@ionsCoverage, $tokens[2]);
			}
			if($ionsCoverage[0]+$ionsCoverage[1]+$ionsCoverage[2]+$ionsCoverage[3] > 2.5){
				for(my $i = 0; $i < 27; $i++){
					$line = <RESULT>;
					chomp $line;
					push(@savedLines, $line);
				}
				foreach $line (@savedLines){
					print $line."\n";
				}
				print "\n";
			}	
		}
	}
}