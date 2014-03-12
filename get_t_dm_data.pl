#!/usr/bin/perl
#
# Calculate dM/dt for a halo using a given HaloHistory output file.

use warnings;
use strict;
use File::Spec;

# Do some rather ugly stuff so that we don't have to place the libraries in some
# @INC directory and can just keep all the script files in one place
my (undef, $modpath, undef) = File::Spec->splitpath($0);
push @INC, File::Spec->canonpath($modpath);
require AMIGA::HaloAnalyzer; import AMIGA::HaloAnalyzer qw(parse_datafile);


#
# Main program code
###################

(@ARGV == 1)
	or die "Arguments: halo index";

my $halo = $ARGV[0];


# Initialize the halos system
my $analyzer = AMIGA::HaloAnalyzer->new(
	data_path=>File::Spec->curdir()
);

# Get the mass evolution data
my %mevo = $analyzer->get_mass_evolution($halo);

# Print it to stdout
foreach my $row (@{$mevo{$halo}}) {
	print join(" ", @$row),"\n";
}
