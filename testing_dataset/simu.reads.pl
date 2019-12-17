#!/usr/bin/perl

#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
#

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <output.prefix>  [num=1e7] [min.size=10] [max.size=200] [read.cycle=100]\n\n";
	exit 2;
}

## set seed, then every time it will give the exactly same result
srand( 7 );

our @nuc = ( 'A', 'C', 'G', 'T' );
our $mutRate = 0.01;
## for testing purpose, using the adapters from the Nextera kits
our $adapter_1 = "CTGTCTCTTATACACATCT";
our $adapter_2 = "AGATGTGTATAAGAGACAG";

my $num   = $ARGV[1] || 1e7;
my $min   = $ARGV[2] || 10;
my $max   = $ARGV[3] || 200;
my $cycle = $ARGV[4] || 100;

## quality is NOT considered for compaison purpose
my $qual  = 'h' x $cycle;

open FQ1, ">$ARGV[0].read1.fq" or die( "$!" );
open FQ2, ">$ARGV[0].read2.fq" or die( "$!" );

foreach my $i ( 1..$num ) {
	my $size = int( rand()*($max-$min) + $min );

	my $read1 = join('', map{ $nuc[int rand @nuc] } (1..$size));
	my $read2 = reverse $read1;
	$read2 =~ tr/ACGT/TGCA/;
	my $m1 = mut($adapter_1);
	my $m2 = mut($adapter_2);
	#print STDERR "$read1+$m1\n$read2+$m2\n";
	$read1 .= $m1;
	$read2 .= $m2;
	if( length($read1) < $cycle ) {
		my $tail = $cycle - length($read1);
		$read1 .= 'A' x $tail;
		$read2 .= 'T' x $tail;
	}

	print FQ1 "\@$i:$size/1\n", substr($read1, 0, $cycle), "\n+\n$qual\n";
	print FQ2 "\@$i:$size/2\n", substr($read2, 0, $cycle), "\n+\n$qual\n";
}

close FQ1;
close FQ2;


sub mut {
	my $raw = shift;
	my $len = length($raw);

	my $mut = '';
	for(my $s=0; $s<$len; ++$s) {
		my $r = substr( $raw, $s, 1 );
		if( rand() < $mutRate ) {	## add a mutation
			while( 1 ) {
				my $m = $nuc[ int rand @nuc ];
				if( $m ne $r ) {
					$mut .= $m;
					last;
				}
			}
		} else {
			$mut .= $r;
		}
	}

	return $mut;
}


