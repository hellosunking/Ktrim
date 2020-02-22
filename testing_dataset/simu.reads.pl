#!/usr/bin/perl

#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
#

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <output.prefix> [num=1e7] [min.size=10] [max.size=200] [read.cycle=100] [error.rate=0.01]\n\n";
	exit 2;
}

## set seed, then every time it will give the exactly same result
srand( 7 );

our @nuc = ( 'A', 'C', 'G', 'T' );
## for testing purpose, using the adapters from the Nextera kits
our $adapter_1 = "CTGTCTCTTATACACATCT";
our $adapter_2 = "AGATGTGTATAAGAGACAG";

my $num   = $ARGV[1] || 1e7;	## number of fragments
my $min   = $ARGV[2] || 10;		## minimum fragment size
my $max   = $ARGV[3] || 200;	## maximum fragment size
my $cycle = $ARGV[4] || 100;	## read cycle
my $eRate = $ARGV[5] || 0.01;	## sequencing error rate

## quality is NOT considered for compaison purpose
my $qual  = 'h' x $cycle;

open FQ1, ">$ARGV[0].read1.fq" or die( "$!" );
open FQ2, ">$ARGV[0].read2.fq" or die( "$!" );

foreach my $i ( 1..$num ) {
	my $size = int( rand()*($max-$min) + $min );

	my $read1 = join('', map{ $nuc[int rand @nuc] } (1..$size));
	my $read2 = reverse $read1;
	$read2 =~ tr/ACGT/TGCA/;
	my $m1 = mut($adapter_1, $eRate);
	my $m2 = mut($adapter_2, $eRate);
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
	my $eRate = shift;

	my $len = length($raw);

	my $mut = '';
	for(my $s=0; $s<$len; ++$s) {
		my $r = substr( $raw, $s, 1 );
		if( rand() < $eRate ) {	## add a mutation
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

