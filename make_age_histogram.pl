#!/usr/bin/perl
#
# Creates a relative amount vs formation time z value histogram for different
# mass bins.


use warnings;
use strict;
use File::Spec;

# Do some rather ugly stuff so that we don't have to place Halos.pm in some
# @INC directory and can just keep all the script files in one place
my (undef, $modpath, undef) = File::Spec->splitpath($0);
push @INC, File::Spec->canonpath($modpath);
require AMIGA::HaloAnalyzer; import AMIGA::HaloAnalyzer qw(parse_datafile);


#
# Main program code
###################

(@ARGV == 2)
	or die "Arguments: processed halodump file, time unit (0: z, 1: gyr)";

my $hdfile = $ARGV[0];
my $use_gyr = $ARGV[1];

my (undef, $datadir, undef) = File::Spec->splitpath($hdfile);
$datadir = File::Spec->curdir() unless ($datadir);

# Initialize the analyzer library. Assume that the data directory is the same
# where the halodump file is.
print "Using $datadir as data path\n";
my $analyzer = AMIGA::HaloAnalyzer->new(
	data_path=>$datadir
);

my @data = &parse_datafile($hdfile);

my @histogram = $analyzer->formation_time_histogram([@data], 
	1.0e8, 1.0e9,
	1.0e9, 1.0e10,
	1.0e10, 1.0e11,
	1.0e11, 1.0e12, 
	1.0e12, 1.0e13,
	1.0e13, 1.0e14,
	1.0e14, 1.0e15,
);

# Add zeroes for those redshift values that had no corresponding halos.
my @zvals = @{$analyzer->{z_values}};
foreach my $bin (@histogram) {
	foreach my $z (@zvals) {
		if (!(exists $bin->[3]{$z})) {
			$bin->[3]{$z} = 0;
		}
	}
}

# Print results in a gnuplot friendly format
for (my $i=0; $i < @histogram; $i++) {
	my $bin = $histogram[$i];
	printf "# Bin %d -- limits (%g, %g) -- halos %d\n",
		$i, $bin->[0], $bin->[1], $bin->[2];
	foreach my $key (sort keys %{$bin->[3]}) {
		my $amount = $bin->[3]{$key};
		print $use_gyr ? $analyzer->z_to_gyr($key) : $key
		      ," ", # time
	       	      $bin->[2] ? $amount/$bin->[2] : 0, "\n"; 
		      # halo fraction
	}
	# Add two extra linebreaks to signify block end for gnuplot
	print "\n\n";
}
