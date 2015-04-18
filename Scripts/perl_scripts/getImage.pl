#!/usr/bin/perl
#grep a list of plots from specplot correspond to the list of IDs
use strict;
my $dir = $ARGV[0];
open(IDS, $ARGV[1]) || die "cannot open file";
my $outdir = $ARGV[2];
my %table;
my @files = `ls $dir`;


while(<IDS>){
    my $id = $_;
    chomp $id;
    $table{$id} = 1;
}

foreach my $file(@files){
    #print "file: $file\n";
    chomp $file;
    $file =~ /^([A-Z0-9\.\+\[\]]+)/;
    my $fileid = $1;
    #print "id: $fileid\n";
    if($table{$fileid}){
	#print $file."\n";
	my $cmd="cp $dir\/$file $outdir";
	print $cmd."\n";
	`$cmd`;
    }
}

