#!/usr/bin/perl
#
# Adds main halo and subhalo lifetime columns to the given 
# halodump_halos_with_subhalos file.

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

(@ARGV == 3)
or die 
<<"HALT";
Arguments: IN OUT CRITERION 
IN: Name of the halodump_halos_with_subhalos.dat file
OUT: Name of the desired output file
CRITERION: A constant defining the birth of a halo.
HALT
;

my ($hdfile, $outfile, $criterion) = @ARGV;

# Guess that the halodata and the HaloHistory binary reside in the same
# directory as the halodump file
my (undef, $dir, undef) = File::Spec->splitpath($hdfile);
$dir = './' unless ($dir);

# Initialize the analyzer
my $analyzer = AMIGA::HaloAnalyzer->new(
	#verbose => 1,
	data_path => $dir, 
	formation_criterion => $criterion);
print "Using formation criterion M(z_form) = ", $criterion, " x M(z=0)\n";


# Get the halodump and modify it for our needs
my @halodump = &parse_datafile($hdfile);
my @modified_dump = $analyzer->modify_halodump(\@halodump);

# Open the output file and print the data there
open (OUT, ">$outfile")
	or die "Couldn't open file for writing: $outfile";

foreach (@modified_dump) {
	print OUT join(" ", @$_)."\n";
}

close OUT;

