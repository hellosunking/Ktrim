#!/usr/bin/perl

#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
#

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <raw.fq> <trimmed.fq> [min.size=36] [seq.cycle=100]\n\n";
	exit 2;
}

my $cutSize  = $ARGV[2] || 36;
my $seqCycle = $ARGV[3] || 100;

my %raw;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	/^\@(\d+):(\d+)\//;
	$raw{$1} = $2;
	<IN>;
	<IN>;
	<IN>;
}
close IN;

my %trim;
open IN, "$ARGV[1]" or die( "$!" );
while( <IN> ) {
	/^\@(\d+):/;
	my $id = $1;
	my $seq = <IN>;
	$trim{$id} = length($seq)-1;
	<IN>;
	<IN>;
}
close IN;

my $perfect  = 0;
my $mistrim  = 0;	## reads that should be trimmed but kept
my $overtrim = 0;	## usually tail-hits
my $overkill = 0;	## the reads are overtrimmed and get discarded
my $error    = 0;	## the reads is trimmed at a wrong position

foreach my $id ( keys %raw ) {
	if( $raw{$id} < $cutSize ) {	## too short reads, should be discarded
		if( exists $trim{$id} ) {
			++ $mistrim;
		} else {
			++ $perfect;
		}
	} else {	## long reads, should be kept
		unless( exists $trim{$id} ) {
			++ $overkill;
			next;
		}
		if( $raw{$id} >= $seqCycle ) {	## very long read
			if( $trim{$id} == $seqCycle ) {
				++ $perfect;
			} else {
				++ $overtrim;
			}
		} else {
			if( $trim{$id} == $raw{$id} ) {
				++ $perfect;
			} elsif( $trim{$id} < $raw{$id} ) {
				++ $overtrim;
			} elsif( $trim{$id} >= $seqCycle-2) {	## seems a tail-trim
				++ $mistrim;
			} else {
				++ $error;
				#print STDERR "$id\n";
			}
		}
	}
}

print "Correct\t$perfect\n",
	  "Over-trim\t",  $overtrim + $overkill, "\n",
	  "Missed-trim\t", $mistrim + $error, "\n";

