#!/usr/bin/perl
#
# Uses a halodump file that has been processed by the add_columns script to
# display data of the number of subhalos that all halos born at a certain
# redshift have.

use warnings;
use strict;
use File::Spec;

# Do some rather ugly stuff so that we don't have to place the libraries in some
# @INC directory and can just keep all the script files in one place
my (undef, $modpath, undef) = File::Spec->splitpath($0);
push @INC, File::Spec->canonpath($modpath);
require AMIGA::HaloAnalyzer; import AMIGA::HaloAnalyzer qw(parse_datafile);


# Main program code
###################

die "Arguments: halodump file, time unit (0 = z, 1 = gyr)" if (@ARGV != 2);
my $hfile = $ARGV[0];
my $use_gyr = $ARGV[1];

# Initialize the halos system
my $analyzer = AMIGA::HaloAnalyzer->new(
	data_path=>File::Spec->curdir()
);

# Constants defining the structure of the file we process.
# Where to find the formation redshift of the main halo
my $zindex = 74;

# Where to find the formation redshift of the subhalo (not needed ATM though).
my $zindex_sub = 105;

# Fetch the data
my @data = &parse_datafile($hfile);

# Extract the number of subhaloes and the formation time for each main halo
my @result;
my $current_halo = -1;
foreach my $row (@data) {
	if ($row->[0] != $current_halo) {
		$current_halo = $row->[0];
		push @result, [
			$use_gyr ? $analyzer->z_to_gyr($row->[$zindex]) : 
				$row->[$zindex],
		       	1];
	}
	else {
		$result[$#result]->[1]++;
	}
}

# Sort
sub by_time {
	$a->[0] <=> $b->[0];
}

@result = sort by_time @result;


# Print what we got (could do something else, too)
foreach my $row (@result) {
	print "$row->[0] $row->[1]\n";
}

