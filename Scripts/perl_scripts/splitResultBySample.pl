#!/usr/bin/perl
open(SAMPLES, $ARGV[0]) || "cannot open sample list file";
my $resultFile = $ARGV[1];

@samples = <SAMPLES>;

for my $sample(@samples){
    chomp $sample;
    my $outFile = $sample;
    $outFile =~ s/\.[a-zA-Z]+$/\.swathout/; 
    $cmd = "cat $resultFile | awk '{if(\$1 ~ /$sample/){print}}' > $outFile";
    print "$cmd\n";
    `$cmd`;
}
