#!/usr/bin/perl

my $header = "#SpectrumFile\tScan#\tAnnotation\tProtein\tCharge\tMQScore\tLength\tTotalPRMScore\tMedianPRMScore\tFractionY\tFractionB\tIntensity\tNTT\tp-value\tF-Score\tDeltaScore\tDeltaScoreOther\tRecordNumber\tDBFilePos\tSpecFilePos\tPrecursorMZ\tPrecursorMZError\n";


open(RESULT, $ARGV[0]) || die "cannot open result file\n";


print $header;
while(<RESULT>){
    my $line = $_;
    if($line !~ /\#/){
	chomp $line;
	my @tokens = split(/\t/, $line);
	print "$tokens[0]\t$tokens[1]\t$tokens[7]\t$tokens[7]\t$tokens[6]\t";
	for(my $i = 0; $i < 5; $i++){
	    print "1\t";
	}
	print "$tokens[11]\t";
	for(my $i = 0; $i < 11; $i++){
	    print "1\t";
	}
	print "\n";
    }
}
