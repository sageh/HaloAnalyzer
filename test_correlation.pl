#!/usr/bin/perl


use warnings;
use strict;

use AMIGA::HaloAnalyzer qw(parse_datafile);
use AMIGA::Vector;
use AMIGA::Correlation qw(ls_correlation);

die "arguments: 3 x N datafile" if (@ARGV != 1);


# Read the data.
my @bf = &parse_datafile($ARGV[0]);

# Vectorize.
for (my $i=0; $i < @bf; $i++) {
	$bf[$i] = AMIGA::Vector->new(@{$bf[$i]});
}

# Calculate correlation
my @result = &ls_correlation(
	#AMIGA::Vector->new(0, 0, 0),
	#AMIGA::Vector->new(1, 1, 1),
	low_bound => AMIGA::Vector->new(-1, -1, -1),
	high_bound => AMIGA::Vector->new(2, 2, 2),
	data => \@bf,
	#[[0.0, 0.25], [0.25, 0.5], [0.5, 0.75], [0.75, 1.0]]);
	bins => [
	[0.01, 0.03162], [0.03162, 0.1], [0.1, 0.3162], [0.3162, 1.0],
	[1.0, 3.162], [3.162, 10.0]
	],
	external_prg => './ls_correlator'
);

print "Result:\n";
foreach my $row (@result) {
	print "$row->[0][0] $row->[0][1] $row->[1] $row->[2]\n";
}

