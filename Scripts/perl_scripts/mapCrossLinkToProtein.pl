#!/usr/bin/perl
#map putative cross-link pepties back to protein to get 
#overall evidence of cross-linked residues

open(RESULT, $ARGV[0]) || die "cannot open MXDB search results";
open(PROTEIN, $ARGV[1]) || die "cannot open protein file";

my $protein="";
my $start=0;
my $pepIndex1=10;
my $pepIndex2=12;
while(<PROTEIN>){
    if($_ =~ />/){
	$start=1;
	next;
    }
    if($start){
	chomp $_;
	$protein=$protein.$_;
    }
    
}
#print "protein is: ".$protein."\n";

while(<RESULT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    my $pep1 = $tokens[$pepIndex1];
    $pep1 =~ /([0-9\.\+]+)/;
    #print "pep1 is: $pep1\n";
    my $offset1 = $1;
    my $position1 = index($pep1, $offset1);     
    $pep1 =~ s/[0-9\.\+\.]+//g;
    #print length($pep1)."\n";
    my $pep2 = $tokens[$pepIndex2];
    $pep2 =~ /([0-9\.\+]+)/;
    my $offset2 = $1;
    my $position2 = index($pep2, $offset2); 
    $pep2 =~ s/[0-9\.\+\.]+//g;
    #print "peptide1: $pep1\n";
    #print "protein: $protein\n";
    my $proteinPos1 = index($protein, $pep1);
    my $proteinPos2 = index($protein, $pep2);  
    my $link1 = $proteinPos1+$position1;
    my $link2 = $proteinPos2+$position2;
    print "ID: ".$tokens[$pepIndex1]."--".$tokens[$pepIndex2]."\t"."map to protein position\t";
    if($proteinPos1 < $proteinPos2){
	print $proteinPos1."\t".$proteinPos2."\t";
    }else{
	print $proteinPos2."\t".$proteinPos1."\t";
    }
    if($link1 < $link2){
	print $link1."\t".$link2."\n";
    }else{
	print $link2."\t".$link1."\n";
    }
    
    
}
