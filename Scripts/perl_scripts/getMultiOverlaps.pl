#!/usr/bin/perl
#compute overlap for multiple list of results
my $numOfRuns = 3;
for(my $f = 0.01; $f < 0.10; $f = $f+0.01){
    my $cmd = "cat search_report.txt | awk '{if(\$1 ~ /PepList/ && \$5  <= ".$f." && \$2 ~ /msplit_[578]/ && \$7 !~ /UPS/) print \$3}' | sort | uniq -c | sort -k1,1nr | awk '{if(\$1 >= $numOfRuns) print}' | wc -l";
    my $count = `$cmd`;
    print "overlap-count at $f $count\n";
}
