#!/usr/bin/perl
use strict;

open(RESULT, $ARGV[0]) || die "cannot open result file";

while(<RESULT>){
    my $line = $_;
    chomp $line;
    if($line =~ /M\[147\]/){
	$line =~ s/\[147\]/\+15\.9945/g;
    }
    if($line =~ /C\[160\]/){
	$line =~ s/\[160\]//g; #no need for explicit PTM takend care in residue mass for C
    }
    if($line =~ /C\[143\]/){ #Pyro-cmC according to library description
	$line =~ s/\[143\]/\-17\.0257/g; #no need for explicit PTM takend care in residue mass for C
    }
    if($line =~ /C\[339\]/){ #ICAT_heavy 
	$line =~ s/\[339\]/\+179\.1357/g; #236.157185-57.021464
    }
    if($line =~ /C\[330\]/){  #ICAT_Light 
	$line =~ s/\[330\]/\+170\.1055/g; #227.126991-57.021464
    }
    if($line =~ /n\[43\](\w)/){
	$line =~ s/n\[43\]/\+42\.010565/g;
    }
    if($line =~ /Q\[111\]/){
	$line =~ s/\[111\]/\-17\.026549/g;
    }
    print $line."\n";
    
}
