#!/usr/bin/perl
#remove wrong annotation (i.e those do not match MS1 precursor mass)
#from MSPLIT annotation
open(IDS, $ARGV[0]);
open(ANNOT, $ARGV[1]);

%idTable=();
@ids = <IDS>;
@idTable{@ids} =1;

while(<ANNOT>){
    my $line = $_;
    chomp $line;
    my @tokens = split(/\s+/, $line);
    if(exists $idTable{$tokens[0]."\n"}){
		#print "found\t fixing annotation\t";
		print $tokens[0]."\t".$tokens[1]."\t".$tokens[1]."\n";
    }else{
		print $line."\n";
    }
}
