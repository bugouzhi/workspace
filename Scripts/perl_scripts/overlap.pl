#!/usr/bin/perl
open($fh1, $ARGV[0]);
open($fh2, $ARGV[1]);
open($fh3, $ARGV[2]);
my $index1=$ARGV[3];
my $index2=$ARGV[4];
my $index3=$ARGV[5];

@array1 = parseFile($fh1, $index1);
@array2 = parseFile($fh2, $index2);
@array3 = parseFile($fh3, $index3);

sub parseFile{
    my $file = $_[0];
    my $index = $_[1];
    @array = ();
    while(<$file>){
	#print $_;
	#chomp $_;
	my @tokens = split(/\t/, $_);
	#print $tokens[0]."\n";
	my $pep = $tokens[$index];
	if($pep=~ /^.\..+\..$/){
	    #print "$pep\n";
	    $pep = substr($pep, 2, length($pep)-4);
	    #print "$pep\n";
	}
	#print $pep;
	push(@array, $pep);
    }
    return @array;
}

#print $#array1."\n";
#print $#array2."\n";
#print $#array3."\n";



%table1 = ();
%table2 = ();
%table3 = ();

@table1{@array1} = 1;
@table2{@array2} = 2;
@table3{@array3} = 3;

%combinetable = ();
@combinetable{@array1}=1;
@combinetable{@array2}=2;
@combinetable{@array3}=3;

%overlap = ();

foreach $item(keys %combinetable){
    #print substr($item, 1, length($item)-2)."\t";
    my $pep = $item;
    $pep =~ s/\n//g;
    print $pep."\t";
    my $code = "";
    if(exists $table1{$item}){
	$code=$code."1 ";
    }else{
        $code=$code."0 ";
    }

    if(exists $table2{$item}){
	$code=$code."1 ";
    }else{
	$code=$code."0 ";
    }
    if(exists $table3{$item}){
	$code=$code."1 ";
    }else{
	$code=$code."0 ";
    }
    if(exists $overlap{$code}){
	$count = $overlap{$code};
	$count=$count+1;
	$overlap{$code} = $count;
    }else{
	$overlap{$code} = 1;
    }
    print "$code\n";
}
#allow some spacing
print "\n";
#print overlap stat summary
print "Overlap stats summary:\n";
foreach my $code(keys %overlap){
    print $code."\t".$overlap{$code}."\n";
}
