#! /usr/bin/perl
use warnings;
use strict;

unless ($#ARGV == 0) {
print "usage: <SAM from STDIN> max_number_of_hits\n";
exit;		
}

my $max = int(shift or die);


## run through each sam line and retain only if number of hits less than specified number

while (my $line = <STDIN>) {
	if ($line =~ /^@/) {
		print STDOUT $line;
	}
	if ($line =~ /NH:i:(\d+)/) { 
		if ($1 <= $max ) {
			print STDOUT $line;
		}
	}
	else {}
}
