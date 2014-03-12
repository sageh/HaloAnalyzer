#!/usr/bin/perl
# 
# Adds a formation time column to a z=0 _halos file.

use warnings;
use strict;
use File::Spec;
use IO::Handle;

# Do some rather ugly stuff so that we don't have to place Halos.pm in some
# @INC directory and can just keep all the script files in one place
my (undef, $modpath, undef) = File::Spec->splitpath($0);
push @INC, File::Spec->canonpath($modpath);
require AMIGA::HaloAnalyzer; import AMIGA::HaloAnalyzer qw(parse_datafile);


# Main program code
###################

(@ARGV == 3)
or die 
<<"HALT";
Arguments: IN OUT CRITERION 
IN: Name of the halos file
OUT: Name of the desired output file
CRITERION: A constant defining the birth of a halo.
HALT
;

my ($hdfile, $outfile, $criterion) = @ARGV;

# Set autoflushing if it is for some reason vital to have the stdout output
# in real time.
#autoflush STDOUT 1;

# Guess that the halodata and the HaloHistory binary reside in the same
# directory as the halodump file
my (undef, $dir, undef) = File::Spec->splitpath($hdfile);
$dir = './' unless ($dir);

# Initialise the analyzer library
my $analyzer = AMIGA::HaloAnalyzer->new(
	data_path=>$dir,
	formation_criterion=>$criterion,
	verbose=>1); # Set to 0 to suppress excess stdout blabber

# Print some general info
print "Reading data from $dir.\n";
print "Using formation criterion M(z_form) = ", $criterion,
" * M(z=0)\n";

# Get the halos_file
my @halodata = &parse_datafile($hdfile);

# DEBUG: Do some testing with a smaller set.
#my @halos = (2000..2002);
#print "Searching for halos: ".join(" ", @halos)."\n";
#my %hh = $analyzer->halohistory(@halos);
#print "Got halohistories for halos:\n";
#foreach my $key (keys %hh) {
#	print "$key ";
#}
#print "\n";
#
#foreach my $halo (@halos) {
#	my @age = $analyzer->find_age($hh{$halo});
#	print "halo: $halo age: $age[1]\n";
#}
#exit;

# Get the halo histories
print "Gathering halo histories...\n";
my %halohistory = $analyzer->halohistory(0..$#halodata);

# Go through the halos file line by line, and add the formation time column
for (my $i=0; $i < @halodata; $i++) {
	# Get formation time for this halo 
	print "Finding formation time for halo $i... ";
	my $age = ($analyzer->find_age($halohistory{$i}))[0];
	print "z = $age\n";
	push @{$halodata[$i]}, $age;
}

# Open the output file and print the data there
open (OUT, ">$outfile")
	or die "Couldn't open file for writing: $outfile";

foreach my $row (@halodata) {
	print OUT join(" ", @$row)."\n";
}

close OUT;

