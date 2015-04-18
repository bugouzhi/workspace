#!/usr/bin/perl


for(my $i = 0; $i <= $#ARGV; $i++){
    open(LIB, $ARGV[$i]) || die "cannnot open $ARGV[$i]\n";
    if($i>0){
	<LIB>;
    }
    while(<LIB>){
	print $_;
    }
}
