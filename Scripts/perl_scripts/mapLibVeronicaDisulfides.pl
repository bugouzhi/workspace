#!/usr/bin/perl
#map whether peptide is a disulfide library peptide from veronica 
use strict;
open(RESULT, $ARGV[0]) || die "cannot open file";
my $pepInd1 = 11;
my $pepInd2 = 13;

while(<RESULT>){
    my $line =$_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    #mapping first peptide
    my $match = 0;
    if($tokens[$pepInd1] =~ /K*[AW]*[DE]*F*[VSHY]A[DY]SC.+VA[KR]$/){
    #if($tokens[$pepInd1] =~ /F[VSHY]A[DY]SC.+VA[KR]/){
	$match = 1;
    }

    if($tokens[$pepInd1] =~ /[TW]*A*[LE]*H*[FV]SC.+VT[PSGY]F[KR]$/){
	$match = 1;
    }

    if($tokens[$pepInd1] =~ /[WA]*V*K*[FL]*C.+[DE]T[VSGY]FA[KR]$/){
	$match = 1;
    }
    
    my $match2 = 0;
    #mapping second peptide
    if($tokens[$pepInd2] =~ /K*[AW]*[DE]*F*[VSHY]A[DY]SC.+VA[KR]$/){
	$match2 = 1;
    }

    if($tokens[$pepInd2] =~ /[TW]*A*[LE]*H*[FV]SC.+VT[PSGY]F[KR]$/){
	$match2 = 1;
    }

    if($tokens[$pepInd2] =~ /[WA]*V*K*[FL]*C.+[DE]T[VSGY]FA[KR]$/){
	$match2 = 1;
    }

    if($match==0){
	$tokens[$pepInd1] = "r".$tokens[$pepInd1];
    }
    
    if($match2==0){
	$tokens[$pepInd2] = "r".$tokens[$pepInd2];
    }
    
    for(my $i=0; $i <= $#tokens; $i++){
	print $tokens[$i]."\t";
    }
    print "\n";
}
