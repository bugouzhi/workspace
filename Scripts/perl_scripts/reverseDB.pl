#!/usr/bin/perl
#reverse a protein DB
open(DB, $ARGV[0]);
my $header="";
my $prevLines = "";

while(<DB>){
	my $line = $_;
	chomp $line;
	if($line =~ />/){
	    if(length($prevLines) > 1){
		print ">X_".substr($header,1)."\n";
		print substr(reverse($prevLines), 1, length($prevLines))."\n";
	    }
	    $header = $line;
	    $prevLines = "";
	}else{
	    $prevLines = $prevLines.$line;
	}
}
print ">X_".substr($header,1)."\n";
print substr(reverse($prevLines), 1, length($prevLines))."\n";
$header = $line;
$prevLines = "";
