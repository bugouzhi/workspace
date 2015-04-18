#!/usr/bin/perl
#compute distance to nearrest tryptic residues for non-tryptic termini

open(FASTA, $ARGV[0]) || die "cannot open fasta files";
my $protein = "";
my $beginProtein=0;
while(<FASTA>){
    my $line = $_;
    chomp $line;
    if($line =~ /^>/){
	$beginProtein=1;
	$protein=$protein."*";
    }elsif($beginProtein){
	$protein=$protein.$line;
    }
}

open(RESULT, $ARGV[1]) || die "cannot open result file";

while(<RESULT>){
    my $line = $_;
    #print "line is $line\n";
    chomp $line;
    my @tokens = split(/\t/, $line);
    my $pep = substr($tokens[7], 2, length($tokens[7])-4);
    $pep =~ s/[0-9\.\+\-]//g;
    my $ind = index($protein, $pep);
    if($ind < 0){
	#print "index is $ind\t$pep\t$tokens[7]\n";
    }
    for(my $i = $ind; $i >= 0; $i--){
	if(substr($protein, $i, 1) =~ /[CSTVKR\*]/ || $i == 0){
	    my $long = substr($protein, $i, length($pep)+($ind-$i));
	    print "peptide\t".$pep."\toffset: ".($i-$ind)."\t".$long."\n";
	    last;
	}
    }
}
