#!/usr/bin/perl
open(ANNON1, $ARGV[0]) || die "cannot open file 1";
open(ANNON2, $ARGV[1]) || die "cannot open file 2";

my %table = ();
while(<ANNON1>){
	my $line = $_;
	chomp $line;
	my @tokens = split(/\t/, $line);
	$table{$tokens[0]} = $tokens[1];
}

while(<ANNON2>){
	my $line = $_;
	chomp $line;
	my @tokens = split(/\t/, $line);
	if(exists $table{$tokens[0]}){
		if($table{$tokens[0]} =~ $tokens[1]){
			print "match\t$tokens[0]\t$table{$tokens[0]}\t$tokens[1]\n";
		}else{
			print "not-match\t$tokens[0]\t$table{$tokens[0]}\t$tokens[1]\n";
		}
	}else{
		print "Not-found\t$tokens[0]\t$tokens[1]\n";
	}
}


