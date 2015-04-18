#!/usr/bin/perl
#generate a peakview library based on peptides overlap
#with spectral ibrary
open(LIBPEPS, $ARGV[0]) || die "cannot open spectral library peptides";
open(PVLIB, $ARGV[1]) || die "cannot open peakview lib";
%table={};
$pepInd = 6;
while(<LIBPEPS>){
    $pep = $_;
    chomp $pep;
    $table{$pep} = 1;
} 

while(<PVLIB>){
    $line = $_;
    @tokens = split(/\t/, $line);
    #print "peptide $tokens[$pepInd]\n";
    if($table{$tokens[$pepInd]}){
	print $line;
    }
}
