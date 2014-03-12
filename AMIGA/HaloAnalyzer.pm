# HaloAnalyzer.pm. Written in 2006 by Pauli Pihajoki
package AMIGA::HaloAnalyzer;

use warnings;
use strict;
use Carp;
use File::Glob ':glob';
use File::Spec;
use File::Temp;

use AMIGA::MergerTree;
use Tools::Vector;
use AMIGA::Correlation qw(calculate_correlation);

our(@ISA, @EXPORT, @EXPORT_OK, $VERSION);

use vars qw($VERSION);
$VERSION = 0.10;

require Exporter;
@ISA = qw(Exporter);
#@EXPORT = qw(); 
@EXPORT_OK = qw(parse_datafile get_mass_interval get_parameter_interval);

#
# Class attributes. These can only be accessed via accessors.
#############################################################

# The instance specific attributes, with given default values. 
my %InstanceAttributes = (
	# Initialisation options
	########################


	# Paths and file names 
	data_path => File::Spec->curdir(),
	binary_path => File::Spec->curdir(),
	halohistory_binary => 'HaloHistory',
	correlator_binary => 'correlator',
	component_correlator_binary => 'compcorr',
	redshift_wildcard => '%z',
	haloindex_wildcard => '%i',
	simufile_template => 'SIMU.z%z.AHF',
	halohistory_output_template => 'output_from_HaloHistory_%i.dat',
	halohistory_input_template => 'halohistory_input_XXXXX',
	simu_logfile => 'SIMU.log',
	simu_dumpfile => 'SIMU.z0.000.AHF.halodump_halos_with_subhalos.dat',

	# File format
	halodump_columns => 74,
	modified_halodump_columns => 140,

	# Physical parameters
	box_size => 10,
	formation_criterion => 0.50,
	hubble_constant => 0.7,

	# Behaviour control
	verbose => 0,

	# Parameters derived at runtime
	###############################
	z_values => undef,
	z_a_t_relation => undef 
);

my %GlobalAttributes = (
	# Column names and of a _halos file
	column_names => [
	"npart", "nvpart", "Xc", "Yc", "Zc", "VXc", "VYc", "VZc", "Mvir",
	"Rvir", "Vmax", "Rmax", "sigV", "lambda", "Lx", "Ly", "Lz", "a", "Eax",
	"Eay", "Eaz", "b", "Ebx", "Eby", "Ebz", "c", "Ecx", "Ecy", "Ecz",
	"ovdens", "Redge", "nbins"],
);

# Generate global accessors
for my $datum (keys %GlobalAttributes) { 
	no strict "refs";    
	*$datum = sub {
		use strict "refs";    
		my ($class, $newvalue) = @_;
		$GlobalAttributes{$datum} = $newvalue if @_ > 1;
		return $GlobalAttributes{$datum};
	}
}

# Generate accessors for instance attributes
for my $datum (keys %InstanceAttributes) { 
	no strict "refs";    
	*$datum = sub {
		use strict "refs";    
		my ($self, $newvalue) = @_;
		$self->{$datum} = $newvalue if @_ > 1;
		return $self->{$datum};
	}
}

#
# Initialization
################ 
sub new {
	my $class = shift;

	# Copy a set of default attribute values
	my $self = {};
	%$self = %InstanceAttributes;

	# Bless ourselves
	bless ($self, $class);

	# Check that the argument list reduces into a hash
	if (@_ % 2) {
		confess "new: arguments: hash table";
	}
	else {
		# It does, so get the arguments, using the accessors, and
		# complain harshly about unknown options
		my %arg = @_;
		foreach my $key (keys %arg) {
			if (exists $self->{$key}) {
				$self->$key($arg{$key});
			}
			else {
				croak "$class: Unknown option $key";
			}
		}
	}

	# Use the templates to figure out what redshift values we have data for,
	# and store the result.
	$self->{z_values} = [$self->find_z_values];

	# Get the redshift - scale parameter - lookback time relation from
	# the SIMU.log file.
	$self->{z_a_t_relation} = [$self->get_time_data];

	# TODO: Make accessors for all new runtime derived attributes, by the
	# same method as for InstanceAttributes

	return $self;
}

###########
# Methods #
###########

#
# Input creators for the AMIGA binaries
#######################################

# Create input the AMIGA MergerTree binary. Returns a string with the correct
# input for the data found in the path.
sub create_mergertree_input {
	my $self = shift;
	my $path = $self->{data_path};
	my $file_template = $self->{simufile_template}."_particles";

	# Construct a filename glob and find all matching _particles files
	my $glb = File::Spec->catfile($path, $file_template);
	$glb =~ s/$self->{redshift_wildcard}/*/g;
	my @partfiles = &bsd_glob($glb);
	if (@partfiles == 0 || !$partfiles[0]) {
		carp "Didn't find any _particles files with template ",
			$file_template, "\n";
		return undef;
	}

	# Construct a regexp for extracting z from the filename
	my $re = $file_template;
	$re =~ s/$self->{redshift_wildcard}/(\\d+\.\\d+)/g;

	# Sort ascending by z (this needs to be a reference to a subroutine,
	# so that the the by_z subroutine and the $re variable don't fall
	# out of sync).
	my $by_z = sub { 
		my ($az, $bz);
		($az) = ($a =~ /$re/);
		($bz) = ($b =~ /$re/);
		$az <=> $bz;
	};
	@partfiles = sort $by_z @partfiles;

	# Generate first those parts of the output that will be
	# identical for each index

	# Number of _particles files
	my $result;
	$result = scalar(@partfiles)."\n";

	# Particles file names
	foreach my $h (@partfiles) {
		$result .= "$h\n";
	}

	# mtree file names
	foreach my $h (@partfiles[0..$#partfiles-1]) {
		$h =~ s/_particles$/_mtree/;
		$result .= "$h\n";
	}

	return $result;
}

# Create input the AMIGA HaloHistory binary.
# Returns an array of strings, each string containing the correct HaloHistory
# input for each halo index given as an argument. The strings are in the
# same order as argument halo indexes.
#
# Arguments: halo indexes
sub create_halohistory_input {
	my $self = shift;
	# Check that the arguments are valid.
	my $path = $self->{data_path};
	my $t = $self->{simufile_template}."_halos";
	my @i = @_;
	my $redshift_wildcard = $self->{redshift_wildcard};
	my $haloindex_wildcard = $self->{redshift_wildcard};

	# Construct a filename glob and find all matching _halos files
	my $glb = File::Spec->catfile($path, $t);
	$glb =~ s/$redshift_wildcard/*/g;
	my @halofiles = &bsd_glob($glb);
	if (@halofiles == 0 || !$halofiles[0]) {
		carp "Didn't find any _particles files with template ",
			$t, "\n";
		return undef;
	}

	# Construct a regexp for extracting z from the filename
	my $re = $t;
	$re =~ s/$redshift_wildcard/(\\d+\.\\d+)/g;
	# Sort ascending by z (this needs to be a reference to a subroutine,
	# so that the the by_z subroutine and the $re variable don't fall
	# out of sync).
	my $by_z = sub { 
		my ($az, $bz);
		($az) = ($a =~ /$re/);
		($bz) = ($b =~ /$re/);
		$az <=> $bz;
	};
	@halofiles = sort $by_z @halofiles;

	# Generate first those parts of the output that will be
	# identical for each index

	# Redshifts
	my $num_redshifts = @halofiles;
	my $redshifts;
	foreach my $h (@halofiles) {
		my ($z) = ($h =~ /$re/);
		$redshifts.= "$z\n";
	}

	# Halo file names
	my $halofiles;
	foreach my $h (@halofiles) {
		$halofiles .= "$h\n";
	}

	# mtree_idx file names
	my $mtreefiles;
	foreach my $h (@halofiles[0..($#halofiles-1)]) {
		$h =~ s/_halos$/_mtree_idx/;
		$mtreefiles .= "$h\n";
	}

	# Loop over each index, generating output
	my @result;
	foreach my $index (@i) {
		my $input;
		my $hh_out = $self->{halohistory_output_template};
		$hh_out =~ s/$haloindex_wildcard/$index/g;

		# First, print the index to track
		$input = $index."\n";

		# Then print the number of redshifts
		$input .= $num_redshifts."\n";

		# Redshifts
		$input .= $redshifts;

		# Halo file names
		$input .= $halofiles;

		# mtree_idx file names
		$input .= $mtreefiles;

		# HaloHistory output file name
		my $out = $self->{halohistory_output_template};
		$out =~ s/%i/$index/g;
		$input .= File::Spec->catfile($path, $out);

		# Output string complete, store it
		push @result, $input;
	}

	return @result;
}

#
# Methods for applying AMIGA binaries to the data
#################################################

# Runs AMIGA HaloHistory program for each given halo index.
# Arguments: halo indexes
sub run_halohistory {
	my $self = shift;
	my $binpath = $self->{binary_path};
	my $datapath = $self->{data_path};
	my $t = $self->{simufile_template};
	my @a = @_;

	# Create HaloHistory input for each halo
	my @hh_input = $self->create_halohistory_input(@a);

	# Loop over the inputs and use a temporary file to give to
	# HaloHistory
	foreach my $input (@hh_input) {
		# Create a temporary file in the data directory
		my $tmpfile = new File::Temp(
			TEMPLATE => $self->{halohistory_input_template},
			DIR => $datapath);

		# Store the data in it
		print $tmpfile $input;

		# Run halohistory
		my $runcmd = File::Spec->catfile($binpath, 
			$self->{halohistory_binary})
			." < $tmpfile 1>".File::Spec->devnull();
		print "Going to run: $runcmd\n" if ($self->verbose);
		system($runcmd) == 0
			or carp "Running HaloHistory failed";
	}
}



#
# Methods for creating different perl structures based on the data
##################################################################

# Finds all the z-values that we have simulation data for.
#
# Arguments: data path
sub find_z_values {
	my $self = shift;
	my $path = $self->{data_path};

	# Make a suitable fileglob by using _halos files, because those are
	# essential to have for every redshift.
	my $glb = File::Spec->catfile($path, 
		$self->{simufile_template}."_halos");
	$glb =~ s/$self->{redshift_wildcard}/*/g;
	
	# Get filenames and check that we actually found matching files
	my @fs = &bsd_glob($glb);
	if (@fs == 0 || !$fs[0]) {
		carp "Couldn't extract redshift information from files ",
	       		"at path: ", $path, "\n";
		return undef;
	}

	# For each of these, extract the redshift
	my @zs;
	foreach (@fs) {
		my @z = &get_floats_from_filename($_);
		if (@z > 1) {
			carp "More than one possible redshift value found in ",
			"filename $_:\n".join("\n",@z)."\n";
			return undef;
		}
		elsif (@z == 0) {
			carp "No valid redshift values found in filename $_\n";
			return undef;
		}	
		push @zs, $z[0];
	}

	# Sort by z in ascending order
	@zs = sort {$a <=> $b} @zs;

	return @zs;
}

# Extract the redshift - scale parameter - time correspondency data from a 
# SIMU.log file
# Arguments: -
sub get_time_data {
	my $self = shift;
	my $path = $self->data_path;
	my $fname = $self->simu_logfile;

	unless (open (IN, "<".File::Spec->catfile($path,$fname))) {
		carp "Couldn't open file for reading: $fname"; 
		return undef;
	}

	# Read until we reach the right section, then extract the data
	my $errmsg = "Couldn't find time data from $fname.";
	while (<IN>) {
		unless (defined($_)) {
			carp $errmsg;
			return undef;
		}
		last if (/^output files:/);
	}

	# Next line has just a separator
	unless (<IN> =~ /^=+\s*$/) {
		carp $errmsg." Didn't find a separator.";
		return undef;
	}

	# Then we have the amount of outputs
	my ($num_redshifts) = (<IN> =~ /^\s*(\d+)\s.*$/);
	if ($num_redshifts < 0) {
		carp $errmsg." Didn't find redshifts.";
		return undef;
	}

	# Next line has just the column names;
	unless (<IN> =~ /^\s*z\s+a\s+t\s*$/) {
		carp $errmsg." Didn't find column names.";
		return undef;
	}

	# Then we have the actual data
	my @result;
	for (my $i=0; $i < $num_redshifts; $i++) {
		my $r = <IN>;
		chomp($r);
		$r =~ s/^\s*//;
		my @row = split (/\s+/, $r);
		if (@row != 3) {
			carp $errmsg." Found ".scalar(@row).
				" rows instead of 3.";
			return undef;
		}
		push @result, [@row]; 
	}

	return @result;
}

# Get all _halos data for a given redshift
#
# Arguments: redshift
sub get_all_halodata {
	my $self = shift;
	my $path = $self->{data_path};
	croak "get_all_halodata: arguments: redshift" if (@_ != 1);
	my $z = shift;

	# Open the correct _halos file
	my $fname = $self->{simufile_template}."_halos";
	$fname =~ s/$self->{redshift_wildcard}/$z/;

	# Get the data in an array
	my @result = &parse_datafile(File::Spec->catfile($path, $fname));

	return @result;
}

# Uses MergerTree mtree_idx files to track halos throughout the simulation in
# the same way as the AMIGA HaloHistory program. Returns a hash table keyed by
# halo indexes containg two dimensional arrays of redshift and corresponding
# index.
#
# Arguments: halo indexes
sub track_halos {
	my $self = shift;
	croak "track_halos: arguments: halo indexes" if (@_ < 1);
	my $path = $self->{data_path};
	my @halos = @_;

	# Get all redshifts
	my @zvals = @{$self->{z_values}};

	# Get data from every _idx file
	my @idx_data;
	foreach (my $i=0; $i < @zvals-1; $i++) {
		my $fname = $self->{simufile_template}."_mtree_idx";
		$fname =~ s/$self->{redshift_wildcard}/$zvals[$i]/;
		push @idx_data, [&parse_datafile(
			File::Spec->catfile($path, $fname))];
	}

	# Track each given halo
	my %result;
	foreach my $oldindex (@halos) {
		# Store the original index separately. $oldindex will be
		# overwritten.
		my $original_index = $oldindex;

		# Add the info for the current (z=0) situation
		my @tmparray;
		push @tmparray, [$zvals[0], $oldindex];

		# Check the idx file for each z value (minus one, because we
		# have one less _idx files than _halos files)
		foreach (my $i=0; $i < @zvals-1; $i++) {
			# Find our halo
			print "Tracking halo $oldindex at z = $zvals[$i]"
			if ($self->verbose);
			my $index_cmp = sub {
				my ($x, $y) = @_;
				return $x->[0] == $y;
			};
			my $idx = &array_find_index_ex($idx_data[$i], 
				$oldindex, $index_cmp);
			if (defined $idx) {
				# Found it. Push the z-value and the
				# corresponding halo index and re-iterate 
				$oldindex = $idx_data[$i]->[$idx]->[1];
				push @tmparray, [$zvals[$i+1], $oldindex];
				print " --> halo $oldindex at ",
				"z = $zvals[$i+1]\n" if ($self->verbose);
			}
			else {
				# We didn't find the halo, so it's useless to
				# continue
				print "\nDidn't find halo $oldindex at ",
				"z=$zvals[$i]\n" if ($self->verbose);
				last;
			}
		}

		$result{$original_index} = [@tmparray];
	}

	return %result;
}

# Function in an identical manner to the AMIGA HaloHistory binary. Output the
# _halos data as a hash for each of the tracked halos for each redshift.
#
# Arguments: Either:
# 	1) Hash with:	use_existing => any false or true value 
# 			halos => [array reference to a list of halo indexes]
#	2) List of halo indexes
# use existing: If this argument is a true value, then use existing AMIGA
# HaloHistory output files when available. Default is true.
sub halohistory {
	my $self = shift;
	my $data_path = $self->data_path;

	my $usage_msg = "halohistory: arguments: [use_existing=>{0,1}] ".
			"halos=>[halo indexes]";

	my %args;
	my $use_existing = 1;
        my @halos;

	# Check the arguments.
	if (@_ % 2) {
		# Definitely not a hash table
		@halos = @_;
	}
	else {
		if (@_ == 0) {
			# No arguments at all.
			croak $usage_msg;
		}

		# So we should have a hash then.
		%args = @_;

		# Make sure we have the keys.
		unless (exists $args{use_existing} 
			&& exists $args{halos}) {
			# We didn't, so guess that it's a table
			# of halo indexes after all.
			@halos = @_;
		}
		else {
			# It was a hash, so just copy the data.
			@halos = @{$args{halos}};
			$use_existing = $args{use_existing};
		}
	}

	# For redshift sorting
	sub by_z {
		$a->[0] <=> $b->[0];
	}

	# Prepare for tracking only those haloes that we didn't find ready made
	# halo history files for.
	my @halos_to_track;

	# Get the halohistory for each halo
	my %result;
	my $infile;
	my @tmparray;
	if ($use_existing) {
		HALO: for my $halo (@halos) {
			# Clear the temporary data storage array
			@tmparray = ();

			print "Checking for pre-existing halo history ",
			      "for halo $halo\n" if ($self->verbose);

			# Check if we have AMIGA HaloHistory output available
			# Get the corresponding filename
			$infile = $self->halohistory_output_template;
			$infile =~ s/$self->{haloindex_wildcard}/$halo/g;
			$infile = File::Spec->catfile($data_path, $infile);

			# Check if it already exists. If it does, we assume
			# that it is valid and use it. Otherwise, we will get
			# the data ourselves.
			if (-e $infile) {
				print "Using existing file for halo $halo\n";
				# Open the file and split into an array
				open (IN, "<$infile") 
					or croak "Couldn't open file for ", 
						 "reading: $infile";
				while (<IN>) {
					chomp;
					# Push an array reference so that we
					# get a two dimensional array as a
					# result
					push @tmparray, [split];
				}
				close IN;

				@tmparray = sort by_z @tmparray;

				$result{$halo} = [@tmparray];
				next HALO;
			}
			else {
				print "\nDidn't find existing history ",
				"for halo $halo\n" if ($self->verbose);
				# Didn't find a ready history, so we will have
				# to create it ourselves.
				push @halos_to_track, $halo;
			}
		}
	}
	else {
		# If we don't want to use any pre-existing data, then we will
		# have to generate it all.
		@halos_to_track = @halos;
	}

	# Now handle all those halos for which we didn't use any pre-existing
	# halo history data for.
	my @z_vals;
	my @halodata;
	my %halotrack;

	# Check if we have any halos to process in the first place
	if (@halos_to_track) {
		# We have, so first fetch all z values and halo data for
		# all z.
		@z_vals = @{$self->{z_values}};
		for my $z (@z_vals) {
			push @halodata, [$self->get_all_halodata($z)];
		}

		# Then track the halos across the redshifts.
		%halotrack = $self->track_halos(@halos_to_track);
	}

	# Then perform a search trough the idx data for each halo, and fetch
	# the halo data, in an identical manner to the AMIGA HaloHistory
	# program.
	for my $halo (@halos_to_track) {
		# Clear the temporary data storage array
		@tmparray = ();

		# There was no existing data, or we didn't want to use it.
		print "Building halo history for halo $halo.\n"
			if ($self->verbose);

		# Get the halos data for each redshift and output the correct
		# row, preceded by the redshift value
		my @ht = @{$halotrack{$halo}};
		for (my $i=0; $i < @ht; $i++) {
			push @tmparray, [$ht[$i]->[0], 
				@{$halodata[$i]->[$ht[$i]->[1]]}];
		}

		# Before adding the array to the hash table, also save it to a
		# file for future use.

		# Get the corresponding filename
		$infile = $self->halohistory_output_template;
		$infile =~ s/$self->{haloindex_wildcard}/$halo/g;
		$infile = File::Spec->catfile($data_path, $infile);

		# Write the data.
		# XXX: At the moment, the data written isn't formatted in
		# fixed width columns as by HaloHistory and thus isn't very
		# human readable. This should maybe be rectified?
		print "Writing $infile.\n" if ($self->verbose);
		open (OUT, ">$infile")
			or croak "Couldn't open file for writing: $infile";
		foreach my $row (@tmparray) {
			print OUT join(" ", @$row), "\n";
		}
		close OUT;

		$result{$halo} = [@tmparray];
	}

	return %result;
}

# Get data for the mass evolution of halos as a function of time. Creates
# a hash of halo indexes pointing to 2D arrays with rows consisting of:
# redshift, time (Gyr), M, dz, dt (Gyr), dM, dM/dz, dM/dt
#
# Arguments: halo indexes
sub get_mass_evolution {
	my $self = shift;

	croak "get_mass_evolution: arguments: halo indexes" if (@_ == 0);
	my @halos = @_;

	# Get the halohistories for these halos
	my %hhs = $self->halohistory(@halos);

	# For each halo, go through the data, and construct the mass evolution
	# data set.
	my %result;
	my @hh;
	my ($t1, $t2, $dt, $dz, $dm);
	foreach my $halo (keys %hhs) {
		$result{$halo} = [];
		@hh = @{$hhs{$halo}};
		for (my $i=0; $i < @hh-1; $i++) {
			# Get the times in gigayears corresponding to the two
			# redshifts.
			$t1 = $self->z_to_gyr($hh[$i+1]->[0]);
			$t2 = $self->z_to_gyr($hh[$i]->[0]);
			$dt = abs($t1-$t2);
			$dz = abs($hh[$i]->[0] - $hh[$i+1]->[0]);
			$dm = $hh[$i]->[9] - $hh[$i+1]->[9];

			# Calculate the values and store the row.
			push @{$result{$halo}}, 
			[$hh[$i]->[0], $t2, $hh[$i]->[9], 
			$dz, $dt, $dm, $dm/$dz, $dm/$dt];
		}
	}

	# Return result
	return %result;
}


# Read a _substructure file and create a hash table with halo index pointing
# to an array of subhalo indexes.
#
# Arguments: redshift
sub get_all_subhalo_data {
	my $self = shift;
	my $path = $self->{data_path};
	croak "get_all_subhalodata: arguments: redshift" if (@_ != 1);
	my $z = shift;

	# Open the correct _halos file
	my $fname = $self->{simufile_template}."_substructure";
	$fname =~ s/$self->{redshift_wildcard}/$z/;

	# Get the data in an array
	my @substr = &parse_datafile(File::Spec->catfile($path, $fname));

	# First row should be ignored, after that we have alternating 
	# rows of: 
	# - halo index and amount of subhalos
	# - subhalo indexes

	# From this, construct a hash table with halo index pointing to
	# an array reference of indexes of subhalos of that main halo.
	# Then the number of subhalos can be obtained from using Perl scalar
	# context on the subhalo index array.
	my %result;
	for (my $i=1; $i < @substr-1; $i+=2) {
		$result{$substr[$i]->[0]} = [@{$substr[$i+1]}];
	}

	return %result;
}

	
#
# Methods for manipulating halo data structures
###############################################

# Add the age and formation time halo data columns to a halodump data
# passed in as an array reference.
#
# Arguments: halodump data as an array reference
sub modify_halodump {
	my $self = shift;
	my $path = $self->data_path;

	croak "modify_halodump: arguments: halodump data as an array reference"
       	if (@_ == 0);
	my $halodump = shift;

	# Print some progress information
	print "Finding halo ages and adding columns...\n" if ($self->verbose);
	my $rowcount = 0;

	# Get halodata and substructure data for all redshifts.
	my @hd;
	my @ss;
	for my $z (@{$self->z_values()}) {
		print "Getting substructure and halodata for z = $z\n"
		if $self->verbose();
		push @hd, [$self->get_all_halodata($z)];
		push @ss, {$self->get_all_subhalo_data($z)};
	}

	# Go through the data line by line, and add the new columns
	my $current_mainhalo = $halodump->[0][0];
	my %hhistory = $self->halohistory($current_mainhalo);
	my @mainhalo_cols = $self->find_age($hhistory{$current_mainhalo});
	my %tmphash = $self->track_halos($current_mainhalo);
	my @mainhalo_track = @{$tmphash{$current_mainhalo}};
	my @result;
	foreach my $row (@$halodump) {
		# Print progress dots
		if ($self->verbose) {
			if ($rowcount++ > 68) {
				$rowcount = 0;
				print "\n";
			}
			else {
				print ".";
			}
		}

		# Check if the main halo changed
		if ($row->[0] != $current_mainhalo) {
			$current_mainhalo = $row->[0];
			%hhistory = $self->halohistory($current_mainhalo);
			@mainhalo_cols = $self->find_age(
				$hhistory{$current_mainhalo});

			# We need the mainhalo track
			%tmphash = $self->track_halos($current_mainhalo);
			@mainhalo_track = @{$tmphash{$current_mainhalo}};
	       }

	       # Get the subhalo column data, and then add the whole deal
	       %hhistory = $self->halohistory($row->[1]);
	       my @subhalo_cols = $self->find_age($hhistory{$row->[1]});

	       # In addition to the halo history data, for each subhalo
	       # we also want to find:
	       # - redshift when it became a subhalo
	       # - redshift when it became a subhalo of it's z=0 mainhalo
	       # - halodata from the latter of the two redshifts
	       
	       # For this, first track the halo
	       %tmphash = $self->track_halos($row->[1]);
	       my @ht = @{$tmphash{$row->[1]}};

	       # Then iterate from the end and find the redshifts
	       my ($zsh1, $zsh2) = (undef, undef);
	       my @zhd;
	       for (my $i=$#ht; $i >= 0; $i--) {
		       # Subhalo id
		       my $id = $ht[$i]->[1];

		       print "Checking subhalo status for subhalo $id ",
		       "at redshift ",$self->z_values()->[$i],"\n"
		       if $self->verbose();

		       # Go over subhalo data, and check if this halo is
		       # the subhalo of some halo, but only if we've not
		       # established that yet
		       if (not defined $zsh1) {
			       for my $halo (keys %{$ss[$i]}) {
				       for my $sid (@{$ss[$i]->{$halo}}) {
					       if ($sid == $id) {
						       $zsh1 = $i;
					       }
				       }
			       }
		       }

		       # If the main halo wasn't formed at this time, then
		       # it can't have any subhalos, so skip to next redshift
		       next if ($i >= @mainhalo_track);

		       # Main halo id
		       my $main_id = $mainhalo_track[$i]->[1];
		       print "Checking if the main halo is $main_id\n"
		       if $self->verbose();

		       # Also specifically check, if the subhalo is now the
		       # subhalo of it's z=0 mainhalo
		       #print "Main id: $main_id\n, keys: ",
		       #join(" ", sort keys (%{$ss[$i]})),"\n";
		       if (not defined $zsh2 and 
			       exists $ss[$i]->{$main_id}) {
			       for my $sid (@{$ss[$i]->{$main_id}}) {
				       if ($sid == $id) {
					       $zsh2 = $i;
					       @zhd = @{$hd[$i]->[$id]};
				       }
			       }
		       }

		       # If we have both, then no need to loop anymore
		       last if (defined $zsh1 and defined $zsh2);
	       }

	       # Check the results
	       if (not defined $zsh1) {
		       # If the first redshift is not defined, then the
		       # other one obviously also isn't. Set to 0.
		       $zsh1 = $zsh2 = 0;
		       @zhd = @{$hd[0]->[$ht[0]]};
		       print "zsh1 not defined\n";
	       }

	       if (not defined $zsh2) {
		       # If the second redhift is not defined, then assume
		       # that it's 0. It must be something for we are only
		       # handling halos that in reality are subhalos!
		       $zsh2 = 0;
		       @zhd = @{$hd[0]->[$ht[0]]};
		       print "zsh2 not defined\n";
	       }

	       # Push the relevant info to @subhalo_cols
	       print "Resulting values -- z1: ",
	       $self->z_values()->[$zsh1], " z2: ",
	       $self->z_values()->[$zsh2], "\n" if $self->verbose();
	       push @subhalo_cols, $self->z_values()->[$zsh1], 
	       		$self->z_values()->[$zsh2], @zhd;

	       push @result, [@$row, @mainhalo_cols, @subhalo_cols];
       } 

	return @result;
}


#
# Methods for deriving further results from the data
####################################################

# Finds the formation time of a halo from HaloHistory output that has
# been parsed into an array. Returns an array containing mass at z=0 and
# the whole row for z = formation split in an array.
#
# Arguments: mass criterion, halohistory output (reference to an array)
sub find_age {
	my $self = shift;

	croak "find_age: arguments: reference to halo history data array" 
	if (@_ != 1);
	my $hh_output = shift;

	# Go through the list until we match the mass criterion
	my $found = 0;
	my @result; 
	# Add the mass at z = 0
	#push @result, $hh_output->[0][9];
	foreach my $z (@$hh_output) {
		if ($z->[9] <= 
			$self->{formation_criterion}*$hh_output->[0][9]) {
			# Found some definite time
			$found = 1;
			push @result, @$z;
			last;
		}
	}

	# If we at any time didn't exceed the criterion, then return
	# information for the greatest z value.
	if ($found == 0) {
		push @result, @{$hh_output->[$#$hh_output]};
	}

	return @result;
}

# Create a formation time histogram for main halos using given mass bins.
#
# Arguments: An array reference to processed halodump data, mass bins
# mass bins: These are given as pairs of lower and upper limits. 
sub formation_time_histogram {
	my $self = shift;

	croak "formation_time_histogram: arguments: reference to a processed ",
	      "halodump data file, mass bins" if (@_ < 3);

	my $data = shift;
	my @bins = @_;

	# Check that we have a valid halodump data that has been complemented
	# with the additional data columns
	unless (@{$data->[0]} == $self->modified_halodump_columns) {
		carp "formation_time_histogram: not valid processed ",
		     "halodump data";
		return undef;
	}

	# Check that we have an even number of bin limiters
	if (@bins % 2) {
		carp "formation_time_histogram: odd number of bin delimiters";
		return undef;
	}

	# The resulting data structure is an array containing an array reference
	# for each bin. This array reference in turn contains:
	# - bin lower bound
	# - bin upper bound
	# - total number of halos in this bin
	# - a hash table keyed by z-value pointing to number of halos that have
	#   this formation time
	my @result;

	# Add the bins
	for (my $i=0; $i < @bins - 1; $i+=2) {
		push @result, [$bins[$i], $bins[$i+1], 0, {}];
	}

	# Go through the data, and group the halos into the bins
	my $current_main_halo = -1;
	my ($m, $z);
	foreach my $row (@$data) {
		# Check if have a new main halo. If so, then add it into all
		# the fitting bins
		if ($row->[0] != $current_main_halo) {
			$current_main_halo = $row->[0];
			$m = $row->[10];
			$z = $row->[75];
			foreach my $bin (@result) {
				# Compare mass to bin limits
				if ($bin->[0] <= $m && $m < $bin->[1]) {
					# Increment number of halos in this bin
					$bin->[2]++;

					# Increment the number of halos in this
					# bin, at this redshift
					$bin->[3]{$z}++;
				}
				#print "Bin: ".join(" ",@$bin)."\n";
			}
		}

		# Other rows are just subhalo data, and are skipped
	}

	return @result;
}

# Using a MergerTree mtree file calculates for every halo the fraction of their
# particles that every other halo shares with it. This can be used for example
# to find what halos have merged.
#
# Arguments: redshift
sub shared_particles {
	my $self = shift;
	croak "shared_particles: arguments: redshift" if (@_ != 1);
	my $path = $self->data_path;
	my $z = shift;

	# Construct a filename template
	my $fname = $self->simufile_template."_mtree";
	$fname =~ s/$self->{redshift_wildcard}/$z/g;

	# Read the file, contents per row:
	# receiving halo, recv. halo n, common n, contributing halo, ctr. halo n
	my @data = &parse_datafile(File::Spec->catfile($path, $fname));

	# First find out the particles that a contributing halo at the earlier
	# redshift has ending up in any halo. For this, make first sure that the
	# data is sorted ascending by the contributing halo. Store the result
	# in a hash, with halo index as a key to the data.
	@data = sort {$a->[3] <=> $b->[3]} @data;
	my $chalo = -1;
	my %ctr_data;
	foreach my $row (@data) {
		# Check for contributing halo change
		if ($row->[3] != $chalo) {
			$chalo = $row->[3];
			$ctr_data{$chalo} = [$row->[4], $row->[2]];
			next;
		}

		# Otherwise, just increment the number of particles in
		# the last entry
		$ctr_data{$chalo}->[1] += $row->[2];
	}

	# Sort by the receiving halo
	@data = sort {$a->[0] <=> $b->[0]} @data;

	# Go through the data again, and collect results in the form of: - a
	# hash table with receiving halo index pointing to an array reference
	# containing:
	#   - number of particles
	#   - hash table with contributing halo index pointing to an array
	#     reference containing:
	#     - number of particles in contributing halo
	#     - number of common particles with receiving halo
	#     - number of common particles with any halo
	#     - fraction: common n / receiving halo n
	#     - fraction: common n / contributing halo n
	#     - fraction: common n / particles of contributing halo common 
	#       with any halo
	#       	  
	#my @result;
	my %result;
	$chalo = -1;
	foreach my $row (@data) {
		# Check for receiving halo change
		if ($row->[0] != $chalo) {
			$chalo = $row->[0];
			$result{$chalo} = [$row->[1], +{}];
			#push @result, [$chalo, $row->[1], []];
		}

		# Then add data for a new contributor halo to the
		# last entry. For this we need the contributor data from the
		# hash we created earlier. 
		my $parts_in_any_halo = $ctr_data{$row->[3]}->[1];
		${$result{$chalo}->[1]}{$row->[3]} = [
			$row->[4], $row->[2], $parts_in_any_halo,
			$row->[2]/$row->[1], $row->[2]/$row->[4],
			$row->[2]/$parts_in_any_halo];
	}

	return %result;
}

# Construct a merger tree (as in Wechsler et al., 2001) with some criteria for 
# the progenitor halo properties. Return the data as a tree like data
# structure for easy top down (z=0 ->) evaluation.
#
# Arguments: inclusion criterion, maxlevels, halo indexes at z=0, 
#
# inclusion criterion: an array reference containing the following data:
#   - minimum amount of particles that must go from halo A to halo B to consider
#     A a progenitor of B
#   - minimum amount of those particles of halo A that end up in halos
#     (as opposed to ending up free) that must end up in halo B to consider
#     A a progenitor of B
#   - minimum mass that a halo must have, in order to include it in the
#     tree at all
# maxlevels: maximum number of redshift levels to traverse back from z=0
sub merger_tree {
	my $self = shift;
	my $path = $self->data_path; 

	croak "merger_tree: arguments: inclusion criterion, maxlevels, halo indexes, " 
		if (@_ < 3);
	my $c = shift;

	croak "merger_tree: arguments: need 3 inclusion criterion parameters" 
		if (@$c != 3);

	# Find the redshifts
	my @z_vals = @{$self->{z_values}};

	# set maxlevels if used
	my $maxlevel = $#z_vals;
	my $ml = shift;
	croak "merger_tree: arguments: maxlevels < 0" if ($ml < 0);
	if ($ml < $maxlevel-1) {
		$maxlevel = $ml;
	}

	# Glomp the rest of the arguments into an array containing the halo indexes
	my @z0_indexes = @_;

	# We will need _halos data and particle sharing information for
	# each redshift.
	my %halodata;
	my %sharedata;
	foreach my $z (@z_vals[0..$maxlevel-1]) {
		$halodata{$z} = [$self->get_all_halodata($z)];
		$sharedata{$z} = {$self->shared_particles($z)};
	}


	# The resulting data structure will be recursive tree structure as
	# follows:
	# - Hash table with halo id at z=0 pointing to an array reference:
	#   - array reference of the halo data at z_0
	#   - hash with progenitor halo ids pointing to an array references:
	#     - array reference of progenitor halo data at z_1
	#     - hash table with 2nd gen. progenitor halo id pointing to an array
	#       reference:
	#       - array reference of 2nd gen. progenitor halo data at z_2
	#	- hash table with 3rd gen. progenitor halod id:
	#		.
	#		.
	#		.
	# 
	# Et cetera. So the basic structure is an array reference containing
	# the halo index, halo data, and an array reference containing the same
	# data for each progenitor halo. In essence a singly linked list.

	# The final resulting data structure will be built in this.
	my %result;

	# Repeat the following procedure for each halo:
	foreach my $index_z0 (@z0_indexes) {

	# Add starting data
	$result{$index_z0} = [$halodata{$z_vals[0]}->[$index_z0], {}];

	# Do a recursive search through the data
	my $z_idx = 0;
	&find_and_add_nodes($result{$index_z0}, $index_z0, $c,
		\@z_vals, \%halodata, \%sharedata, $z_idx, $ml);

	# Recursively search for halos matching criterion
	sub find_and_add_nodes {
		# Found progenitor halos go here
		my $ph = shift;
		# Find progenitors for this halo index
		my $hi = shift;
		# Using this criterion
		my $c = shift;
		# z values
		my @z_vals = @{shift;};
		# Halo data for each z
		my $hd = shift;
		# Particle sharing data for each z
		my $sd = shift;
		# Index to z_vals table for this level
		my $zi = shift;
		# Maximum z level
		my $mz = shift;

		# Stop if we are at the last possible redshift
		if ($zi >= $mz-1) {
			return;
		}

		# Check if we have progenitor candidates for this halo, and
		# get them, if we do. Stop otherwise.
		if (exists $$sd{$z_vals[$zi]}{$hi}) {
			my %pc = %{$$sd{$z_vals[$zi]}{$hi}[1]};

			# Go through each progenitor candidate, and check
			# for match against the criterion.
			foreach my $key (keys %pc) {
				# Get the halo data for the contributing
				# halo
				my @chalo = @{$$hd{$z_vals[$zi+1]}[$key]};
				
				# Compare the relative particle amounts and
				# the mass
				if ($pc{$key}->[4] >= $c->[0] &&
					$pc{$key}->[5] >= $c->[1] &&
					$chalo[8] >= $c->[2]) {
					# Halo passed the progenitor criteria.
					# Add it to the hash.
					$$ph[1]{$key} = [\@chalo, {}];
				}
			}
		}
		else {
			return;
		}


		foreach my $key (keys (%{$$ph[1]})) {
			&find_and_add_nodes($$ph[1]{$key}, $key, $c,
				\@z_vals, $hd, $sd, $zi + 1, $mz);
		}
	}

	# Then make it a MergerTree proper by invoking the constructor
	$result{$index_z0} = new AMIGA::MergerTree($result{$index_z0});
	}

	# Now we should have a hashful of MergerTrees in the %result hash,
	# ready for use.
	return %result;
}

# Calculate the value of Lanzy-Szalay (1993) correlation function for the given
# halodata, using the given scale bins. If transverse scale bins are given, then
# assume that both parallel and transverse correlations are to be calculated.
# Returns a two dimensional array reference containing the bins and the
# corresponding correlation function value.
#
# Arguments: hash with:
# halodata => halodata (arrayref)
# bins => (parallel) scale bins (arrayref)
# trbins => transverse scale bins (arrayref) (optional)
# field => precalculated binomial field (optional)
sub correlation {
	my $self = shift;

	croak "correlation: argument didn't reduce to a hash table" 
	if (@_ % 2);
	my %args = @_;

	my $halos = $args{halodata};
	my $bins = $args{bins};

	my $trbins = undef;
	if (exists $args{trbins}) {
		$trbins = $args{trbins};
	}

	# If a precalculated binomial field was specified, use it.
	my $binfield = undef;
	if (exists $args{field}) {
		$binfield = $args{field};
	}

	if (defined $trbins) {
		carp "correlation: No parallel scale bins specified." 
			unless(@$bins);
		carp "correlation: No transverse scale bins specified." 
			unless(@$trbins);
	}
	else {
		carp "correlation: No scale bins specified." unless(@$bins);
	}
	carp "correlation: No halodata." unless(@$halos);

	# The calculations will take a long time, so we can read the data as
	# we go, instead of reading all halodata at once, to save some memory.
	# In addition, the correlation algorithms want vectorized data and
	# don't care about all the other data fields, so we'll have to cut
	# off the extra.
	my @halodata = @$halos;
	my @inputdata;

	# Our data is distributed in 3 dimensions, between (0, 0, 0) and
	# (b, b, b) where b is the box size.
	my $low_bound = Tools::Vector->new(0, 0, 0);
	my $high_bound = Tools::Vector->new(
		$self->box_size, $self->box_size, $self->box_size); 

	# Prepare an array reference where we store the data
	my @result;

	# Fetch data for the halos specified, and vectorize it
	for (@halodata) {
		push @inputdata, Tools::Vector->new(@{$_}[2..4]);
	}
	print "inputdata: ",scalar(@inputdata),"\n"
	if ($self->verbose());

	# Calculate the correlation and store
	# it in the array.
	# XXX: Use the fast external correlator
	my %corr_args = (
		low_bound => $low_bound, 
		high_bound => $high_bound, 
		data => \@inputdata, 
		bins => $bins
	);
	if (defined $binfield) {
		$corr_args{field} = $binfield;
	}
	if (defined $trbins) {
		$corr_args{trbins} = $trbins;
		$corr_args{external_prg} = File::Spec->catfile(
			$self->binary_path(),
			$self->component_correlator_binary());
	}
	else {
		$corr_args{external_prg} = File::Spec->catfile(
			$self->binary_path(),
			$self->correlator_binary());
	}

	@result = &calculate_correlation(%corr_args);

	return @result;
}

#
# Methods for miscellaneous calculations
########################################

# Convert from redshift to gigayears 
# (not the lookback time: z = 0 -> t = age of the universe).
# Arguments: redshift
sub z_to_gyr {
	my $self = shift;
	croak "z_to_gyr: arguments: redshift" if ((@_) != 1);
	my $z = shift;

	# We could also calculate the look-back time ourselves from the
	# definition of T = (1/H0)int_0^z [(1+z)**2*sqrt(z\Omega_0 + 1)]^-1 dz
	# but this is probably too cumbersome.

	# Get the SIMU.log (z, a, t) relation array
	my @time_data = @{$self->{z_a_t_relation}};

	# Go through the data and try to match the given redshift by finding
	# the value closest to the given one.

	## Index of the candidate
	#my $cnd = 0;
	#for (my $i=0; $i < @time_data; $i++) {
	#	if (abs($time_data[$i]->[0] - $z) 
	#		< abs($time_data[$cnd]->[0] - $z)) {
	#		# Found a better match
	#		$cnd = $i;
	#	}
	#}
	#
	## Before returning the value, check that the difference is not greater
	## than some epsilon. For now, use 0.05.
	#my $epsilon = 0.05;
	#
	#if (abs($time_data[$cnd]->[0] - $z) > $epsilon) {
	#	# Didn't find it. Carp and return undef to signify error;
	#	carp "Couldn't find redshift $z from $self->simu_logfile";
	#	return undef;
	#}
	#else {
	#	# Return the value scaled by the hubble constant to get
	#	# gigayears in standard units.
	#	return $time_data[$cnd]->[2]/$self->hubble_constant;
	#}

	# XXX: The previous algorithm has problems with simulation restarts,
	# because these can cause outputs at very closely separated redshifts.
	# This in turn causes two different redshifts to have the same
	# value in gigayears, which isn't good. So, we use a linear
	# interpolation from closest points instead.

	# Linearly interpolate between closest values found.
	# Note: The @time_data table is sorted by descending redshift value.
	my $li= 0;
	for (my $i=0; $i < @time_data-1; $i++) {
		if ($z < $time_data[$i]->[0]) {
			$li = $i;
		}
	}

	my ($lx, $hx) = ($time_data[$li]->[0], $time_data[$li+1]->[0]);
	my ($ly, $hy) = ($time_data[$li]->[2], $time_data[$li+1]->[2]);
	my $t = ($hy-$ly)/($hx-$lx)*($z-$lx) + $ly;

	return $t/$self->hubble_constant();
}

# Convert the realspace coordinates of halodata to redshift coordinates
# by taking into account the peculiar velocity of each halo for some vantage
# point.
#
# Arguments: halodata (arrayref), observation point (3 components, arrayref)
sub realspace_to_z_space {
	my $self = shift;

	croak "realspace_to_z_space: arguments: halodata, obervation point"
		if (@_ != 2);
	my $halodata = shift;
	my $point = shift;

	# Go through the halodata and calculate results. Use Tools::Vectors
	# for convenience.
	my @result;
	$point = Tools::Vector->new(@$point);
	for (@$halodata) {
		# Get the velocity and relative position vectors
		my $vel = Tools::Vector->new(@$_[5..7]);
		my $pos = Tools::Vector->new(@$_[2..4]);
		#print "pos: $pos vel: $vel\n";
		my $pos_rel = $pos - $point;

		# We also need the normalized relative position vector
		my $pos_n = $pos_rel*(1.0/$pos_rel->length());
		#print "pos_rel: $pos_rel \npos_n: $pos_n\n";

		# Peculiar velocity v_pec = dot(r_hat,v)*r_hat. So, shift the
		# position by v_pec/H_0 = v_pec/(h * 100 km/s)
		$vel = ( $pos_n->dot($vel) / 
				($self->hubble_constant()*100.0) ) * $pos_n;

		#print "Shifting with: $vel\n";
		$pos = $pos + $vel;
		#print "Shifted pos: $pos\n";

		# Store result
		my @row = @$_;
		@row[2..4] = $pos->to_array();
		push @result, [@row];
	}

	return @result;
}


######################
# Helper subroutines #
######################

# Extracts all unsigned decimal form floating point numbers found in a file name
# Argument: file name
sub get_floats_from_filename {
	my $fname = shift;
	my @retval;
	
	while ($fname =~ /(\d+\.\d+)/g) {
		push @retval, $1;
	}

	return @retval;
}

# Reads in a file containing two dimensional whitespace separated data,
# and splits it into a two dimensional array
#
# Arguments: filename (complete with path)

# XXX: Use explicit memoizing of arguments to trade RAM for speed. For
# modify_halodump at least this results in a huge speedup.
my %memoize;
sub parse_datafile {
	croak "parse_datafile: arguments: filename (with path)" if (@_ != 1);
	my $infile = shift;

	# If we have already read this file once, just return what we
	# parsed the last time.
	if (exists $memoize{$infile}) {
		return $memoize{$infile};
	}

	# Otherwise, parse the file.

	# Open the file and split into an array
	open (IN, "<$infile") 
		or croak "Couldn't open file for reading: $infile";
	my @result;
	while (<IN>) {
		chomp;
		# Ignore comments
		next if (/^#/o);

		# Push an array reference so that we get a two dimensional
		# array as a result
		push @result, [split];
	}
	close IN;

	# Memoize it for future use.
	$memoize{$infile} = \@result;

	return @result;
}

# Searches for the first entry in an array equal to argument. Equality is
# defined by a user given comparator. Return index if found, undef if not
# found.
sub array_find_index_ex {
	my $arrayref = shift;
	my $value = shift;
	my $cmpfunc = shift;

	for (my $i=0; $i < @$arrayref; $i++) {
		return $i if (&$cmpfunc($arrayref->[$i], $value));
	}

	# Found none
	return undef;
}

# Extracts a the portion of a halo file that has halos in the given
# mass interval. Returns the portion as an arrayref, as well as the indexes
# corresponding to the highest and lowest mass halo included. 
#
# Arguments: halodata (arrayref), mass low bound, mass high bound
sub get_mass_interval {
	croak "get_mass_interval: arguments: halodata, mass high bound, ".
		"mass low bound" if (@_ != 3);

	return &get_parameter_interval($_[0], 8, $_[1], $_[2]);
	
	# XXX: Use get_parameter_interval internally
	#my ($halodata, $lb, $hb) = @_;
	#if ($lb >= $hb) {
	#	carp "get_mass_interval: low bound is greater than or equal to"
	#	." high bound";
	#	my $tmp = $lb;
	#	$lb = $hb;
	#	$hb = $tmp;
	#}
	#
	## DEBUG:
	##printf "halodata: %d lb: %e hb: %e\n", scalar(@$halodata), $lb, $hb;
	#
	## Do a binary search for both high and low bounds, using the fact that
	## the halodata is in descending mass order.
	#my $li = $#$halodata/2;
	#my $hi = $li;
	#
	#my $l = $#$halodata;
	#my $h = 0;
	#while ($l-$h > 1) {
	#	#print "h $h l $l hi $hi\n";
	#	if ($halodata->[$hi][8] < $hb) {
	#		$l = $hi;
	#		$hi = int(($h + $hi)/2);
	#	}
	#	elsif ($halodata->[$hi][8] > $hb) {
	#		$h = $hi;
	#		$hi = int(($l + $hi)/2);
	#	}
	#}
	## Add one to high index, because interger division will round towards
	## lower index.
	#$hi++;
	#
	#$l = $#$halodata;
	#$h = 0;
	#while ($l-$h > 1) {
	#	#print "h $h l $l li $hi\n";
	#	if ($halodata->[$li][8] < $lb) {
	#		$l = $li;
	#		$li = int(($h + $li)/2);
	#	}
	#	elsif ($halodata->[$li][8] > $lb) {
	#		$h = $li;
	#		$li = int(($l + $li)/2);
	#	}
	#}
	#
	## DEBUG:
	##print "Found li $li hi $hi\n";
	##print "Masses $halodata->[$li][8] $halodata->[$hi][8]\n";
	#
	#return ([@$halodata[$hi..$li]], $hi, $li);
}

# Get those halos from the given halo data with the given parameter in the
# given range.
#
# Arguments: halodata (arrayref) parameter index, parameter low bound,
# parameter high bound.
sub get_parameter_interval {
	croak "get_parameter_interval: arguments: halodata, parameter index, ",
	"parameter low bound, parameter high bound" if (@_ != 4);

	my $halodata = shift;
	my $ix = shift;
	my $lb = shift;
	my $hb = shift;

	if ($lb >= $hb) {
		carp "get_parameter_interval: low bound is greater than ",
		       "or equal to high bound";
		my $tmp = $lb;
		$lb = $hb;
		$hb = $tmp;
	}

	# Sort the halodata descending by the chosen parameter
	@$halodata = sort {$b->[$ix] <=> $a->[$ix]} @$halodata;

	# Do a binary search for both high and low bounds, using the fact that
	# the halodata is now in descending order.
	my $li = $#$halodata/2;
	my $hi = $li;
	
	my $l = $#$halodata;
	my $h = 0;
	while ($l-$h > 1) {
		#print "h $h l $l hi $hi\n";
		if ($halodata->[$hi][$ix] < $hb) {
			$l = $hi;
			$hi = int(($h + $hi)/2);
		}
		elsif ($halodata->[$hi][$ix] > $hb) {
			$h = $hi;
			$hi = int(($l + $hi)/2);
		}
	}
	# Add one to high index, because interger division will round towards
	# lower index.
	$hi++;
	
	$l = $#$halodata;
	$h = 0;
	while ($l-$h > 1) {
		#print "h $h l $l li $hi\n";
		if ($halodata->[$li][$ix] < $lb) {
			$l = $li;
			$li = int(($h + $li)/2);
		}
		elsif ($halodata->[$li][$ix] > $lb) {
			$h = $li;
			$li = int(($l + $li)/2);
		}
	}
	
	# DEBUG:
	#print "Found li $li hi $hi\n";
	#print "Values $halodata->[$li][$ix] $halodata->[$hi][$ix]\n";

	return ([@$halodata[$hi..$li]], $hi, $li);;
}


1;

# Use END tag to begin POD documentation

__END__

=head1 NAME

AMIGA::HaloAnalyzer - AMIGA data analysis package

=head1 SYNOPSIS

	use AMIGA::HaloAnalyzer;

	# Or if we want to use the non-OO functions.
	use AMIGA::HaloAnalyzer qw(parse_datafile 
				   get_mass_interval
				   get_parameter_interval);

	# Initialisation
	my $analyzer = AMIGA::HaloAnalyzer->new(
		# Set paths for external files
		data_path => "path/to/data/", 
		binary_path => "path/to/HaloHistory/and/other/binaries", 
		# Set the physical parameters of the data (default 
		# values shown)
		box_size => 10,
		formation_criterion => 0.50,
		hubble_constant => 0.7,
		# Print additional information about executed processes
		verbose => 1);

	# Create a string suitable for giving as input to the AMIGA 
	# MergerTree program, and use it to run MergerTree.
	my $input = $analyzer->create_mergertree_input();
	system("MergerTree < $input > output_file");

	# Create strings suitable for inputting to the AMIGA 
	# HaloHistory binary. Generating one string per given 
	# halo index.
	my @input = $analyzer->create_halohistory_input(@halo_indexes);

	# Directly run AMIGA HaloHistory program for the given halos.
	$analyzer->run_halohistory(@halo_indexes);

	# Get the contents of a _halos file for a given redshift 
	# as a perl 2D array. The redshift values are best obtained 
	# from the internal z_values array.
	my @halodata = $analyzer->get_all_halodata(
		$analyzer->{z_values}->[0]);
	print "Data for redshift = ",$analyzer->{z_values}->[0],"\n";
	for (my $i=0; $i < @halodata; $i++) {
		print "Mass of halo no. $i is $halodata->[$i][8]\n";
	}

	# Tracing the changes in halo index number across 
	# redshifts starting from z = 0.
	my %halotrace = $analyzer->trace_halos(@halo_indexes);
	foreach my $halo (keys %halotrace) {
		print "Backtrace of halo $halo across redshifts:\n";
		foreach my $row (@{$halotrace{$halo}}) {
			print "Redshift $row->[0], index $row->[1]\n";
		}
	}

	# Retrieving halo history data for given halos.
	my %halohistory = $analyzer->halohistory(@halo_indexes);
	foreach my $halo (keys %halohistory) {
		print "Halo history of halo $halo by ",
		"ascending redshift:\n";
		foreach my $row (@{$halohistory{$halo}}) {
			print join(" ", @$row), "\n";
		}
	}

	# Get data for the mass evolution of halos as a function 
	# of time.
	my %mass_evolution = $analyzer->get_mass_evolution(
		@halo_indexes);
	foreach my $halo (keys %mass_evolution) {
		print "Mass evolution of halo $halo\n";
		print "Fields:\n";
		print "redshift, time (Gyr), M, dz, ",
		"dt (Gyr), dM, dM/dz, dM/dt";
		foreach my $row (@{$mass_evolution{$halo}}) {
			print join(" ", @$row), "\n";
		}
	}

	# Parse a halodump_with_subhalos file and add the 
	# formation time column,and halo data columns corresponding 
	# to formation time to it.
	my @halodump = &parse_datafile(
		"halodump_with_subhalos_file");
	my @modified_halodump = $analyzer->modify_halodump(\@halodump);

	# Find out the formation time of a halo by its halohistory,
	# and return the formation redshift and the halo data for
	# that redshift as an array.
	my %halohistory = $analyzer->halohistory(@halo_indexes);
	foreach my $halo (%halohistory) {
		my @data = $analyzer->find_age($halohistory{$halo});
		print "Halo no. $halo formed at redshift $data[0]\n";
	}

	# Retrieve a mass binned formation time histogram for all
	# the mainhalos in given halodump data.
	my @halodump = &parse_datafile("halodump_file");
	my @modified_dump = $analyzer->modify_halodump(\@halodump);
	my @histogram = formation_time_histogram(\@halodump, 
		[1.0e9, 1.0e10], [1.0e10, 1.0e11], [1.0e11, 1.0e12]);
	foreach my $bin (@histogram) {
		print "Bin with lower bound $bin->[0] and ",
		"upper bound $bin->[1] has $bin->[2] halos.\n";
		print "These are divided between redshifts as ", 
		"follows:\n";
		foreach my $redshift (keys %{$bin->[3]}) {
			print "$bin->[3]->{$redshift} halos ",
			"formed at redshift",
			$redshift, "\n";
		}
	}

	# For some redshift z_i, find out how the halos at z_i and
	# z_i+1 share their particles, by using the _mtree files.
	my %sharing = $analyzer->shared_particles(
		$analyzer->{z_values}->[0]);
	foreach my $halo (keys %sharing) {
		print "Halo $halo at redshift ", 
		$analyzer->{z_values}->[0],
		"has $sharing{$halo}->[0] particles in total.\n";
		print "When compared to halos at redshift ",
		$analyzer->{z_values}->[1], " the common particles are ",
		"divided as follows:\n";
		foreach my $contributor (keys %{$sharing{$halo}->[1]}) {
			print join(" ", 
				@{$sharing{$halo}->[1]->{$contributor}}), 
				"\n";
		}
	}

	# Construct a merger tree structure for each given halo, by
	# applying the given progenitor criterion to the particle
	# sharing data.
	# The criterion is: 
	# - minimum fraction of common particles between progenitor
	#   and offspring, 
	# - minimum fraction of progenitor particles
	#   common with halos that must also be common with the
	#   offspring, 
	# - minimum progenitor mass
	my %merger_trees = $analyzer->merger_tree([0.5, 0.7, 1.0e10], 
		@halo_indexes);
	foreach my $tree (keys %merger_trees) {
		print "We have constructed a merger tree ",
		"for halo $tree\n";
		# Convert the tree from single linked (one way) to double
		# linked (two-way), so that one can traverse into both
		# directions from a node.
		$merger_trees{$tree}->double_link();
		# Get a by-level accessor to the tree
		my @view_by_level = $merger_trees{$tree}->level_view();
	}

=head1 DESCRIPTION

C<AMIGA::HaloAnalyzer> is a package for performing various data analysis tasks
on a set of AMIGA halo finder data. These include constructing merger trees
for individual halos, calculating formation times etc.

The HaloAnalyzer package uses the _halos, _mtree and _mtree_idx files from
the data directory to do it's work. The AMIGA MergerTree program should thus
have been run on the data beforehand to generate these files.

=head1 INITIALISATION

=over 4

=item new ( OPTIONS )

This is the constructor for a new AMIGA::HaloAnalyzer object. C<OPTIONS> is
a hash containing the configuration options as key-value pairs. Valid
keys with default values indicated are:

=over 4

=item data_path

The path where the halo finder data files (_halos, _mtree,
_mtree_idx, etc.) can be found.

Default: Current directory (where the calling script was ran from).

=item binary_path

The path where the library searches for useful binaries,
such as HaloHistory and external correlator programs.

Default: Current directory (where the calling script was ran from).

=item halohistory_binary

The name of the AMIGA HaloHistory program.

Default: C<HaloHistory>

=item correlator_binary

Name of the external two-point correlator program.

Default: C<correlator>

=item component_correlator_binary

Name of the external parallel/transverse two-point correlator program.

Default: C<compcorr>

=item redshift_wildcard

Wildcard used in filenames for redshift values.

Default: C<%z>

=item haloindex_wildcard

Wildcard used in filenames for halo indexes.

Default: C<%i>

=item simufile_template

Template for data file names (_halos, etc.)

Default: C<SIMU.z%z.AHF>

=item halohistory_output_template

Template for halo history files.

Default: C<output_from_HaloHistory_%i.dat>

=item halohistory_input_template

Template for input files for the AMIGA
HaloHistory program. The template should conform with the requirements
for C<mktemp> C standard library function. See L<mktemp(3)> or 
L<File::Temp> for reference.

Default: C<halohistory_input_XXXXX>

=item simu_logfile

Name of the AMIGA simulation log.

Default: C<SIMU.log>

=item simu_dumpfile

Name of the halodump_halos_with_subhalos file. From now
on this file will be referred to as a I<halodump file>.

Default: C<SIMU.z0.000.AHF.halodump_halos_with_subhalos.dat>

=item halodump_columns

Number of columns in an unaltered halodump file.

Default: 74

=item modified_halodump_columns

Number of columns in a halodump file in which
formation time column and the corresponding data columns have been added.

Default: 140

=item box_size

Size of the simulation box used for the data set.

Default: 10

=item formation_criterion

The fraction of the C<z = 0> mass of halo that is
considered to be its initial mass. I.e. when tracking a halo backwards
across the redshifts, its formation time is calculated as being at the redshift
when the mass of the halo falls below this parameter.

Default: 0.5

=item hubble_constant

Value of the normalized Hubble constant 
C<h = H/100 km/s/pc>.

Default: 0.7

=item verbose

Controls the verbosity level of the module. Currently only two values are
supported. Value 1 causes the module to output some information of the ongoing
processes to standard output, while value 0 suppresses most output.

Default: 0 

=back

=back

=head2 ACCESSORS

After creating an instance of the C<HaloAnalyzer> class, the configuration
options can be retrieved and changed using accessors. Each option has an
accessor of the same name that sets the option to the first argument and
returns the value set, or returns the current value if no arguments are
given.

	# As default, the instance is created as non-verbose
	my $analyzer = AMIGA::HaloAnalyzer->new();
	print $analyzer->verbose(), "\n"; # Prints 0
	# Set new value.
	$analyzer->verbose(1);
	print $analyzer->verbose(), "\n"; # Prints 1
	# Set new value and return the set value.
	print $analyzer->verbose(0), "\n"; # Prints 0

=head1 ATTRIBUTES

Every configuration option listed at L<new> is also an attribute
of the C<AMIGA::HaloAnalyzer> object. In addition, there are several other
attributes with useful data that are derived runtime, after the object has
been created. All of these attributes have accessors as detailed in 
L</ACCESSORS>.

=over 4

=item z_values

An array of all the redshift values that could be extracted from the files
at L<data_path> using L<simufile_template> as reference. The module obtains
it by calling the L<find_z_values|/find_z_values ( )> method. The values are
stored in
ascending order of numerical value.

Example content: C<["0.000", "0.100", "0.149"]>

=item z_a_t_relation

A 2D array of the correspondence between the redshift, scale parameter and
time values (Gyr/h) of the simulation output times. These are extracted
from the L<simu_logfile> and are used by the L<z_to_gyr|/z_to_gyr ( REDSHIFT )> 
method.

Example content:

	[[7.000, 0.125, 0.554],
	[4.000, 0.200, 1.118],
	[2.250, 0.308, 2.114],
	[1.250, 0.444, 3.584],
	[0.700, 0.588, 5.238],
	[0.400, 0.714, 6.684],
	[0.200, 0.833, 7.996],
	[0.000, 1.000, 9.706]]

=back

=head1 METHODS 

These are methods of an object, and called as C<$object->method()>.

=head2 create_mergertree_input ( )

Creates a string suitable for piping to the AMIGA MergerTree binary. The
function looks at the files residing in L<data_path> and by comparing them
to the L<simufile_template> creates an input string correct for this dataset.
The filenames are stored in the string I<with relative path>.

Arguments: -

Return value: A string value containing MergerTree input.

Example of usage:

	# Construct input
	my $input = $analyzer->create_merger_tree_input();

	# Use it to run MergerTree program
	system("MergerTree < $input > output_file");

=head2 create_halohistory_input ( HALOS )

Creates strings suitable for piping to AMIGA HaloHistory binary. The function
takes an array of halo indexes as an argument and constructs a corresponding
array of input strings.

Arguments: C<HALOS> -- array of halo indexes.

Return value: An array of HaloHistory input strings.

Example of usage:

	# Construct input for some halos
	my @input = $analyzer->create_halohistory_input(0..100);

	# Use these to run HaloHistory program
	for (my $i=0; $i < @input; $i++) {
		system("HaloHistory < $input[$i] > output_number_$i");
	}

=head2 run_halohistory ( HALOS )

Runs the AMIGA HaloHistory binary for each halo index in C<HALOS> and stores
the output according to L<halohistory_output_template>.

Arguments: C<HALOS> -- array of halo indexes.

Return value: -

Example of usage:

	# Run HaloHistory program for some halos, and check that
	# the results got written
	my @halos = (0..10);
	$analyzer->run_halohistory(@halos);
	foreach my $halo (@halos) {
		my $filename = 
			$analyzer->halohistory_output_template();
		my $wildcard = $analyzer->haloindex_wildcard();
		$filename =~ s/$wildcard/$halo/g;
		if (!(-e $filename)) {
			print "File $filename wasn't created!\n";
		}
	}

=head2 find_z_values ( )

Looks at the filenames in the dataset, and retrieves all the redshift values
found by comparing with L<simufile_template>.

Arguments: -

Return value: Array of redshift values.

Example of usage:

	# Fetch the redshift values from simulation 
	# data filenames
	my @z_vals = $analyzer->find_z_values();
	print "Found redshifts: ".join(" ", @z_vals)."\n";

	# However, these are more efficiently obtained from
	# the object itself.
	my @z_vals_2 = @{$analyzer->z_values()}; # These two are
	my @z_vals_3 = @{$analyzer->{z_values}}; # the same thing!
	# Or fetch just the reference
	my $z_vals_ref = $analyzer->z_values();

	# Print smallest value (sorted in ascending order)
	print $z_vals_2[0], "\n";
	print $z_vals_ref->[0], "\n"; # By using the reference

=head2 get_time_data ( )

Goes through L<simu_logfile> and tries to find the table relating the 
simulation output redshifts to time and scale parameter values. 

Arguments: -

Return value: 2D array with redshift, scale parameter and time columns.

Example of usage:

	# Find the values from simulation logfile
	my @data = $analyzer->get_time_data();

	# But they are already fetched when the object was created
	my @data2 = $analyzer->z_a_t_relation(); # more efficient!
	# Print the values
	print "redshift , scale parameter, time (Gyr/h)\n";
	foreach my $row (@data2) {
		print join(" ", @$row), "\n";
	}

=head2 get_all_halodata ( REDSHIFT )

Retrieves the data from a _halos file for the given redshift and parses it
into a 2D array.

Arguments: C<REDSHIFT> -- a string value.

Return value: 2D array with the columns corresponding to a _halos file.

Example of usage:

	# Fetch the data for the smallest redshift (should be 0.000)
	my $redshift = $analyzer->z_values()->[0]; # From reference
	my @halodata = $analyzer->get_all_halodata($redshift);

	# Print out just the number of particles for each halo
	foreach my $row (@halodata) {
		print $row->[0], "\n";
	}

=head2 track_halos ( HALOS )

Uses MergerTree _mtree_idx files to track each halo in C<HALOS> backwards
through time. This is done in a similar manner as in the AMIGA HaloHistory
binary.

Arguments: C<HALOS> -- array of halo indexes.

Return value: A hash table with halo index keys pointing to 2D arrays. The
arrays contain two columns: redshift and the corresponding halo index.

Example of usage:

	# Track every halo that we have at z = 0.
	my @halodata = $analyzer->get_all_halodata(
				$analyzer->z_values()->[0]);
	# One row for each halo, so amount of rows == amount of halos
	my $n = @halodata;
	# Track every halo. Indexes run from 0, so need to subtract one
	my %track = $analyzer->track_halos(0..($n-1));

	# Print out the results
	foreach my $halo (keys %track) {
		print "Track of halo no. $halo\n";
		print "redshift index\n";
		foreach my $row (@{$track{$key}}) {
			print join(" ", @$row), "\n";
		}
	}

=head2 halohistory ( HALO_HASH | HALOS )

Tracks given haloes backwards through time, and retrieves their parameters
at each simulation output redshift. The output is then saved to the directory
L<data_path> by using the template L<halohistory_output_template>. 
The argument can either be either a hash
table, or an array of halo indexes. If the C<HALO_HASH> contains a key
B<use_existing> that points to a false value, then all output files are
rewritten, otherwise existing files will be used.

Arguments: C<HALOS> -- array of halo indexes, I<or> C<HALO_HASH> -- a hash
table with keys B<use_existing> pointing to any true or false value, and
B<halos> pointing to an array reference of halo indexes.

Return value: A hash table with halo indexes pointing to 2D arrays of halo
history data.

Example of usage:

	# Get halohistory data for some halos
	my %hh = $analyzer->halohistory(2, 4, 5, 7);
	# Or by overwriting existing files
	my %hh2 = $analyzer->halohistory(
			use_existing => 0,
			halos => [2, 4, 5, 7]); 

	# Get halo history data for all halos at z = 0, and overwrite old
	# halo history files. NOTE: This takes a relatively long time and
	# is memory intensive! (order of several hundred MB for ~5000
	# halos)
	my $num_halos = scalar($analyzer->get_all_halodata(
				$analyzer->z_values()->[0]));
	my %hh = $analyzer->halohistory(
			use_existing => 0,
			halos => [0..($num_halos-1)]);

=head2 get_mass_evolution ( HALOS )

Calculates data for the mass evolution of each halo in C<HALOS> 
as a function of time.

Arguments: C<HALOS> -- an array of halo indexes.

Return value: A hash table with halo indexes pointing to 2D arrays with
columns: redshift, time (Gyr), M, dz, dt (Gyr), dM, dM/dz, dM/dt

Example of usage:

	# Get mass evolution of the 10 most massive halos, 
	# and print them
	my %data = $analyzer->get_mass_evolution(0..9);
	foreach my $halo (sort keys %data) {
		print "Mass evolution of halo no. $halo\n";
		print "z, t, M, dz, dt, dM, dM/dz, dM/dt\n";
		foreach my $row (@{$data{$halo}}) {
			print join(" ", @$row), "\n";
		}
	}

=head2 get_all_subhalo_data ( REDSHIFT )

Reads _substructure files for the C<REDSHIFT> and returns a hash table
with halo index keys pointing to arrays of subhalo indexes.

Arguments: C<REDSHIFT> -- a string value.

Return value: Hash table with halo indexes pointing to array references
containing indexes of that halos subhalos.

Example of usage:

	# Fetch the data for the smallest redshift (should be 0.000)
	my %shdata = $analyzer->get_all_halodata(
				$analyzer->z_values()->[0]);

	# Print out just the number of subhalos for each halo, and make
	# sure that the halos are listed in ascending order
	foreach my $halo (sort {$a<=>$b} keys %shdata) {
		print "$halo ",scalar(@{$shdata{$halo}}), "\n";
	}

=head2 modify_halodump ( HALODUMP )

Appends additional columns to
a 2D array I<reference> C<HALODUMP> containing the data from 
a L<simu_dumpfile>. The added columns in order are: 

=over 4

=item *

Main halo formation redshift and halo data for that redshift.

=item *

Subhalo formation redshift and halo data for that redshift.

=item *

Redshift when the subhalo first became a subhalo.

=item *

Redshift when the subhalo first became a subhalo of it's C<z=0> main halo and
halo data for that redshift.

=back

Arguments: C<HALODUMP> -- reference to a 2D array containing halodump data.

Return value: 2D array containing the data from C<HALODUMP> appended with
formation redshift and corresponding _halos data columns.

Example of usage:

	# Get the halodump
	my @halodump = &parse_datafile($analyzer->simu_dumpfile());
	# Add columns
	@halodump = $analyzer->modify_halodump([@halodump]);

=head2 find_age ( HALOHISTORY )

Uses halo history data C<HALOHISTORY> and L<formation_criterion> parameter to
find out the formation redshift and the corresponding _halos data.

Arguments: C<HALOHISTORY> -- reference to a 2D array containing output from
the HaloHistory program or C<halohistory> method.

Return value: An array containing the formation redshift and the
corresponding _halos data row.

Example of usage:

	# Fetch halohistory for some halo, and find it's formation time
	# and parameters
	my $halo_id = 1337;
	my %hh = $analyzer->halohistory($halo_id);
	my @age_and_params = $analyzer->find_age($hh{$halo_id});
	print "Halo $halo_id formed at redshift $age_and_params[0]\n";

=head2 formation_time_histogram ( MOD_HALODUMP [, BIN [, BIN, [...]]] )

Calculates formation time histogram data for the given halodump data in
C<MOD_HALODUMP> that must have been processed by the C<modify_halodump>
method beforehand. The main halos are binned in the given mass bins, and
a separate histogram dataset is calculated for each bin.

Arguments:

=over 4

=item *

C<MOD_HALODUMP> -- reference to a 2D array containing halodump data
modified with the C<modify_halodump> method.

=item *

C<BIN> -- array references to a pair of value specifying the bin lower
and upper bounds.

=back

Return value: An array containing array references for each bin. These point
to arrays containing: 

=over 4

=item *

bin lower bound

=item *

bin upper bound

=item *

total number of halos in this bin

=item *

a hash table with redshift keys pointing to how many halos formed at that redshift

=back

Example of usage:

	my @dump = &parse_datafile($analyzer->simu_dumpfile());
	@dump = $analyzer->modify_halodump([@dump]);
	my @hist = $analyzer->formation_time_histogram(\@dump,
		[1.0e10, 5.0e10], [5.0e10, 1.0e11]);

	foreach my $bin (@hist) {
		print "There are $bin->[2] main halos in ",
		"mass range [$bin->[0], $bin->[1]]\n";
	}

=head2 shared_particles ( REDSHIFT )

Uses _mtree files to calculate various results about how the halos at
C<REDSHIFT> and one output time earlier shared particles with each other.

Arguments: C<REDSHIFT> -- a string value. Best obtained from the L<z_values>
array.

Return value: A hash table keyed with halo indexes at C<REDSHIFT> pointing to
a array references, which contain:

=over 4

=item *

Number of particles of the halo (C<N>).

=item *

Hash table with the indexes of the 
I<contributing halos at the previous redshift>
pointing to array references containing:

=over 4

=item *

Number of particles in the contributing halo (C<N_c>).

=item *

Number of shared particles between the two halos (C<N_s>).

=item *

Total number of its particles that the contributing halo shares with I<any>
halo at C<REDSHIFT> (C<N_a>).

=item *

Fraction: C<N_s/N>.

=item *

Fraction: C<N_s/N_c>.

=item *

Fraction: C<N_s/N_a>.

=back

=back

Example of usage:

	# Find out how many potential progenitors a halo at z=0 has.
	my $redshift = $analyzer->z_values()->[0];
	my %sp = $analyzer->shared_particles($redshift);
	foreach my $halo (keys %sp) {
		print "Halo $halo received particles from ",
			scalar(keys %{$sp{$halo}->[1]}), " halos.\n";
	}

=head2 merger_tree ( CRITERION, HALOS )

Constructs a merger tree for each halo in C<HALOS>. This is done by applying
the progenitor criterion C<CRITERION> to the data from the
L<shared_particles|\shared_particles ( REDSHIFT )> method to find the
progenitors of each halo at the earlier simulation output redshift. Then the
same process is performed again recursively on each progenitor until we have a
a complete merger tree for each halo.

Arguments:

=over 4

=item *

C<CRITERION> -- an array reference with three parameters relating to progenitor
candidate A (at an earlier redshift) to B, the halo currently examined 
(at a later redshift):

=over 4

=item *

Minimum fraction of particles in A that must be shared between A and B to
consider A a progenitor of B. I.e. the minimum fraction of particles of A that
must end up in B.

=item *

Minimum amount of those particles in A that end up in I<any> halo at the B's
redshift that must end up in B to consider A a progenitor of B.

=item *

Minimum mass that A must have in order to qualify as a progenitor at all.

=back

All three of these conditions must be met for a halo to be added as
a progenitor.

=item *

C<HALOS> -- an array of halo indexes at redshift z = 0.

=back

Return value: C<AMIGA::MergerTree> object. See L<AMIGA::MergerTree> for
reference.

Example of usage:

	# Construct a merger tree for some halos, with similar 
	# criteria for the particle fractions (0.5, 0.7) as in 
	# Wechsler et al, 2001
	my @halo_indexes = (0, 1, 2);
	my %mt = $analyzer->merger_tree([0.5, 0.7, 2.2e10],
			@halo_indexes);
	# Get a by level view for each halo and print out how
	# many progenitors there are at each z
	foreach my $halo (keys %mt) {
		print "Halo $halo:\n";
		my @lview = $mt{$halo}->level_view();
		for (my $i=0; $i < @lview; $i++) {
			printf "At z=%f there are %d progenitors\n",
			$analyzer->z_values()->[$i],
			scalar(@{$lview[$i]});
		}
	}

=head2 correlation ( ARGS )

Calculate values of the Landy-Szalay correlation function estimator for the
given halo data using the given mass bins. Additionally a precalculated 
binomial random field may be specified. This is recommended because there is
no guarantee of the suitability of the Perl random number generator for
Monte Carlo type simulations.

The arguments should be passed as a hash table C<ARGS> with the following
keys:

=over 4

=item halodata 

Reference to an array of _halos data. Easily obtained by using
L<get_all_halodata|/get_all_halodata ( REDSHIFT )> method for example.

=item bins

Reference to a 2D array with scale low and high bound columns
for each row. If L<trbins> are specified, then parallel bins should be
specified here.

=item trbins

Reference to a 2D array with transverse scale low and high bound columns
for each row. If this key is specified, then the correlation is calculated
separately for parallel and transverse components, taking the observer to
be at point (0, 0, 0) corner of the box. The parallel scale bins are then
specified in L<bins>.

=item field

I<Optional.> Reference to a 2D array with x, y and z columns containing
a binomial random field.

=back

Return value: 2D array containing correlation data. See L<AMIGA::Correlation>
for reference.

Example of usage:

	# Calculate a correlation using an precalculated binomial field
	my @c = $analyzer->correlation(
		halodata => \@halodata,
		bins => [[0.01, 0.1], [0.1, 10], [10, 100]],
		field => \@binomial_field);
	# Print results
	print "r_min r_max ksi error\n";
	for (@c) {
		print join(" ", @$_), "\n";
	}

=head2 z_to_gyr ( REDSHIFT )

Convert redshift to gigayears using the value of L<hubble_constant> attribute,
and the redshift - scale parameter - time correspondency extracted from the
L<simu_logfile>.  Note: Returns time counted from the beginning of the
universe, I<not> the lookback time.

Arguments: C<REDSHIFT> -- a string value. Best obtained from the 
L<z_values> array.

Return value: A scalar value of the time elapsed corresponding to C<REDSHIFT>.
Returns undef in case of error.

Example of usage:

	# First simulation redshift is usually z = 0, and converting
	# this to time should yield the age of the universe.
	my $age = $analyzer->z_to_gyr($analyzer->z_values()->[0]);
	print "Age of the universe is: $age\n";

=head2 realspace_to_z_space ( HALODATA, POINT )

Adjust the halo positions in C<HALODATA> based on the halo velocities and
the observation point given in C<POINT>. The result is the redshift space
distribution of the halos as seen from C<POINT>.

Arguments:

=over 4

=item *

C<HALODATA> -- an array reference to _halos file data.

=item *

C<POINT> -- an array reference containing the x, y and z components of the
observation point

=back

Return value: Copy of the C<HALODATA> structure with the positions
adjusted according to the velocities and the observation point.

Example of usage:

	# Adjust the halo positions to what we would observe in redshift
	# space from (0, 0, 0) corner of the simulation box.
	my @halodata = $analyzer->get_all_halodata(
				$analyzer->z_values()->[0]);
	my @zdata = $analyzer->realspace_to_z_space(\@halodata, [0, 0, 0]);

=head1 UTILITY FUNCTIONS

These functions are static and don't need to be invoked with an object
or classname. However, each utility function must be separately imported
using the standard way: C<use AMIGA::HaloAnalyzer qw(func1 func2 foofunc);>

=head2 parse_datafile ( FILE )

Assumes that C<FILE> contains columns of whitespace separated data. The file
is read in and parsed into a 2D perl array, discarding rows
beginning with #, which is taken as a comment delimiter. Note: No consistency
checking of any kind is done so for example the number of columns could
vary from line to line.

Arguments: C<FILE> -- a string value containing file name complete with 
absolute or relative path.

Return value: 2D Perl array containing the data read in from C<FILE>.

Example of usage:

	my $filename = $analyzer->data_path()."/data.txt";
	my @data = &parse_datafile($filename);

=head2 get_mass_interval ( HALODATA, LOW, HIGH )

Extracts the part of _halos data given in C<HALODATA> that falls between
limiting masses C<LOW> and C<HIGH>. Note: The content in C<HALODATA> is
assumed to be already sorted by mass in descending order, and is not resorted
before searching. The search method used is a binary search, and will thus
fail, if the data is not ordered.

Arguments:

=over 4

=item *

C<HALODATA> -- reference to a 2D array containing _halos data. Easily obtained
by using for example the C<get_all_halodata()> method.

=item *

C<LOW>, C<HIGH> -- scalar values specifying the low and high mass limits.

=back

Return value: An array reference containing:

=over 4

=item *

Splice of the data in C<HALODATA> for halos with mass between C<LOW> and
C<HIGH>.

=item *

Array index corresponding to the highest mass (and thus lowest index value)
in C<HALODATA> that was included in the splice.

=item *

Array index corresponding to the lowest mass (highest index value)
in C<HALODATA> that was included in the splice.

=back

The return value can be expressed in Perl as (actual code):

	return ([@$halodata[$hi..$li], $hi, $li);

Example of usage:

	my @halodata = $analyzer->get_all_halodata(
				$analyzer->z_values()->[0]);
	# Get only halos from certain mass range
	my @slice = &get_mass_interval(\@halodata, 1.0e9, 5.0e9);
	# Replace original data.
	@halodata = @{$slice[0]};

=head2 get_parameter_interval ( HALODATA, INDEX, LOW, HIGH )

Extracts the part of _halos data given in C<HALODATA> that has the data
at column C<INDEX> (starting from 0) fall between
limiting values C<LOW> and C<HIGH>. The data in C<HALODATA> will be sorted
by the that column first, and then the the limiting indexes are found by
binary searching.

Arguments:

=over 4

=item *

C<HALODATA> -- reference to a 2D array containing _halos data. Easily obtained
by using for example the C<get_all_halodata()> method.

=item *

C<INDEX> -- column index of the _halos data. Indexing starts from 0.

=item *

C<LOW>, C<HIGH> -- scalar values specifying the low and high limits.

=back

Return value: The same as with the C<get_mass_interval()> method, except that
the halo data splice will be sorted by the values in the given column, in
descending order.

Example of usage:

	my @halodata = $analyzer->get_all_halodata(
				$analyzer->z_values()->[0]);
	# Get only halos from certain particle amount range
	my @slice = &get_parameter_interval(\@halodata, 
				0, 5000, 15000);
	# Replace original data.
	@halodata = @{$slice[0]};

=head1 NOTES

=head2 REDSHIFTS

The redshifts in L<z_values> array are extracted from the filenames found from
L<data_path>. Thus they are represented as strings so that when needed, they
can be converted back to filenames by simple substitution.

What this means is that passing explicit numeric values for methods that
require redshift values is not a good idea, and most of the time will not
work. The redshifts should always be fetched from the L<z_values> array.

=head2 DATA REPRESENTATION

The difference between methods requiring halo data (found in _halos files), 
halodump data (found in the L<simu_dumpfile>) and modified halodump data
(halodump data processed with the 
L<modify_halodump|/modify_halodump ( HALODUMP )> method) should be 
carefully observed.

All of these data types are manipulated internally as simple Perl 2D arrays,
and the methods don't really do very through validation of data arguments.
Thus passing wrong type of data could succeed without any errors, but the
results will be unpredictably bizarre.

In addition, when reading in data from a file, it is recommended that
the L<parse_datafile|/parse_datafile ( FILE )> function be used. It 
automatically strips comments
from the file. It doesn't check that every row has an equal amount of columns
though, so if extra care is needed, then data consistency should be verified
manually. Giving data with variable column amount as parameter to a method
will at best result in an error or warning, but could go unnoticed and
cause in unpredictable results.

=head1 AUTHOR

Written in 2006 for the benefit of Tuorla Observatory Cosmological Simulations
research group by Pauli Pihajoki.

=head1 COPYRIGHT

Copyright (c) 2006 Pauli Pihajoki. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms
as Perl itself.

=cut

