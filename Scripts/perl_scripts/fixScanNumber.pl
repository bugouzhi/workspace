#!/usr/bin/perl
#inspect seems to not parsing the scan number in the title
#we needed this table to do the mapping
open(SPECTRUMFILE, $ARGV[0]) || die "cannot open mgf file";
open(INSPECTRESULT, $ARGV[1]) || die "cannot open inspect result file";
my $counter = 0;
my @scanNumbers = ();
while(<SPECTRUMFILE>){
	my $line = $_;
	chomp $line;
	if($line =~ /TITLE=Scan Number: (\d+)/){
		my $scanNumber = $1;
		#print "storing $scanNumber\n";
		push(@scanNumbers, $scanNumber);
	}
} 

my $line = <INSPECTRESULT>; #skip header line
while(<INSPECTRESULT>){
	my $line = $_;
	#print "line is $line\n";
	$line =~ /\t(\d+)\t/;
	my $num = $1;
	#print "count is $num\n";
	my $scanNum = $scanNumbers[$num];
	#print "scan #: $scanNum\n";
	$line =~ s/\t\d+\t/\t$scanNum\t/;
	print $line;
}

