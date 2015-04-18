#!/usr/bin/perl
open(PROTEIN, $ARGV[0]);
open(PEPTIDES, $ARGV[1]);
my $protein=""; #one-line proteins
my $proteinName="";
my $peptideCount=0;
my @peptides = ();

while(<PEPTIDES>){
    my $line = $_;
    my @tokens = split(/\t/, $line);
    my $pepitde = $tokens[7];
    $peptide = substr($pepitde, 2, length($peptide)-4);
    $peptide =~ s/[0-9\.\+\-]//g;
    push(@peptides, $peptide);
}
print "Done loading peptides\n";

while(<PROTEIN>){
	my $line = $_;
	chomp $line;
	if($_ =~ /^>/){
	    foreach my $peptide(@peptides){
		if($peptide =~ $protein){
		    $peptideCount++;
		}
	    }
	    print "$proteinName:\t$peptideCount\n";
	    $proteinName = $line;
	    $peptideCount=0;
	}else{
		$protein = $protein.$line;
	}
}


