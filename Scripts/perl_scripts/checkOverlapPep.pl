#!/usr/bin/perl
#check overlapping peptides check whether peptides share suffix prefix

open(PEP1, $ARGV[0]) || die "cannot open peptides file 1";
open(PEP2, $ARGV[1]) || die "cannot open peptides file 2";

@peps = ();
@peps2 = ();

while(<PEP1>){
    $pep = $_;
    chomp $pep;
    push(@peps, $pep);
}

while(<PEP2>){
    $pep2 = $_;
    chomp $pep2;
    push(@peps2, $pep2);
}

print "list size: $#peps $#peps2\n";

$subLen=6;
foreach $pep(@peps){
    foreach $pep2(@peps2){
	#$pref = substr($pep, 0, $subLen);
	#$suf = substr($pep, length($pep)-$subLen, $subLen);
	$pref = $pep;
	$suf = $pep;
	#print "pref $pref suff $suf\n";
	if(length($pep) > 0){
	    #print "here\n";
	    #print "pep2 is $pep2\n";
	    if(($pep2 =~ $pref || $pep2 =~ $suf)){
		print $pep." and ".$pep2." overlap in sequence\n";
	    }
	}
    }
}
