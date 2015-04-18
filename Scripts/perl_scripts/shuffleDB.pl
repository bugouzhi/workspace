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
		my @seqArray = split(//, $prevLines);
		for(my $i = 0; $i < $#seqArray; $i++){
		    my $ind = int(rand($#seqArray));
		    my $temp = $seqArray[$i];
		    $seqArray[$i] = $seqArray[$ind];
		    $seqArray[$ind] = $temp;
		}
		for(my $i = 0; $i < $#seqArray; $i++){
		    print $seqArray[$i];
		}
		print "\n";
		#substr(shuffle($prevLines), 1, length($prevLines))."\n";
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
