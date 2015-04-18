#!/usr/bin/perl
#merge all the result files to generate a single result report
#only for inspect for now
use strict;
my $pattern = $ARGV[0];
my $command = "ls ".$ARGV[0];
print "command is $command\n";
my @files = `$command`;
print "we have $#files files\n";
foreach my $file(@files){	
	chomp $file;
	my $index = -1;
	if($file =~ /yeast_digets(\d+)_/){
		$index = $1;
	}
	
	if($file =~ /MSPLIT_result_yeast(\d+)/){
		$index = $1;
	}
	
    open(RESULT, $file);
    while(<RESULT>){
		my $line=$_;
		my @tokens = split(/\s+/, $line);
		print "yeast".$index."_$tokens[4]"."\t".$line;
    }
} 
