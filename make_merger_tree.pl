#!/usr/bin/perl
#
# Prints out a merger tree in METAPOST language for the given halo.

use warnings;
use strict;
use File::Spec;

# Do some rather ugly stuff so that we don't have to place HaloAnalyzer.pm in
# some @INC directory and can just keep all the script files in one place
my (undef, $modpath, undef) = File::Spec->splitpath($0);
push @INC, File::Spec->canonpath($modpath);

# Now we can import the libraries
require AMIGA::HaloAnalyzer; import AMIGA::HaloAnalyzer;


# Main program code
###################

unless (@ARGV >= 4 && $ARGV[0] < @{AMIGA::HaloAnalyzer->column_names()}) {
	print "Arguments: parameter index (circles), minimum halo mass, ",
	"number of z lookback levels, halo indexes\n",
	"parameter index should be one of the following:\n";
	for (my $i=0; $i < @{AMIGA::HaloAnalyzer->column_names()}; $i++) {
		print "$i => ".AMIGA::HaloAnalyzer->column_names()->[$i]."\n";
	}
	exit 1;
}

my ($p_index, $min_mass, $maxzlevel, @halos) = @ARGV;

# Presume that the script is run from the data directory.
my $datadir = File::Spec->curdir();
$datadir = './' if (!$datadir);

# Initialize the library 
my $analyzer = AMIGA::HaloAnalyzer->new(
	data_path=>$datadir,
);

# Get z values
my @z_vals = @{$analyzer->z_values()};

# Create the merger tree structure
my %tree = $analyzer->merger_tree([0.5, 0.7, $min_mass], $maxzlevel, @halos);

# Fetch all subhalo data
my @shdata;

foreach my $z (@z_vals) {
	push @shdata, {$analyzer->get_all_subhalo_data($z)};
}

#
# Start outputting the picture to stdout
########################################

# Some graphics related constants (how to scale and so on)
my $x_scaling = 0.7;
my $y_scaling = 0.9;
my $circle_scaling = 0.7;

# The minimum distance between nodes on one row.
my $min_separation = 0.1;

# Print any possible preamble
&output_preamble;

# Then process the tree of each halo and draw it
foreach my $halo (@halos) {

# We will want to draw some parameter dependent things. Store the z=0
# parameters for convenience. 
#my @z0_params = @{$tree{$halo}->root_node->[0]};
# XXX: Store the value of the desired parameter at z0 only
my $z0_param = $tree{$halo}->root_node->[0][$p_index];

# Make the tree doubly linked.
$tree{$halo}->double_link();

# Then, create a new viewport to the tree, indexed by the level
my @level_view = $tree{$halo}->level_view();

# Set the width coordinate for the initial halo
push @{$tree{$halo}->root_node}, 0;

# Then go through the tree level by level, and sort firstly by offspring 
# width coordinate, and secondly by the halo mass
my $by_width_and_mass = sub {
	if ($a->[2][3] != $b->[2][3]) {
		return $a->[2][3] <=> $b->[2][3];
	}
	else {
		return $a->[0][8] <=> $b->[0][8];
	}
};

# While calculating the width coordinates, also find out the maximum width so
# that the redshift labels can be placed aesthetically
my $maximum_width = 0;
foreach my $level (@level_view[1..$#level_view]) {
	# Sort first.
	@$level = sort $by_width_and_mass @$level;

	# Total number of halos on this level
	my $total = @$level;

	# We also need the sum of every virial radius
	my $radiussum = 0;
	foreach my $halo (@$level) {
		$radiussum += $halo->[0][9];
	}

	# Then calculate a width coordinate using the combined virial radius as
	# a weight factor and then scaling with some clever function.
	# Also sum up the total width so far.
	my $widthsum = 0;
	push @{$level->[0]}, 0.0;
	for (my $i=1; $i < @$level; $i++) {
		# The constant multiplier here determines the separation between
		# the halo points in units of the (unscaled) parameter circle
		# diameter. For now, the separation is equal to one half the
		# combined parameter circle diameter of the neighbouring halos.
		#
		# So, later when we are drawing the actual parameter circles,
		# we will most probably want to scale them down just a bit,
		# so that the circles will not touch each other.
		my $separation = 
		0.5*($level->[$i-1][0][$p_index] + $level->[$i][0][$p_index])/
		$z0_param + $min_separation;

		push @{$level->[$i]}, $widthsum + $separation;
		$widthsum += $separation;
	}

	# Then normalize the width to some smaller number. For now, normalize
	# to units of sum of separations on a level / virial radius at z=0.
	foreach my $halo (@$level) {
		if ($widthsum != 0) {
			#$halo->[3] *= 1.0/$z0_params[9];
			#$halo->[3] -= $widthsum/(2.0*$z0_params[9]);
			$halo->[3] -= $widthsum/2.0;

			# Check if this coordinate is the greatest so far
			if (abs($halo->[3]) > $maximum_width) {
				$maximum_width = abs($halo->[3]);
			}
		}
		else {
			$halo->[3] = 0; 
		}
	}
}


# Now find out the number of subhalos of the most massive halo on each level.
# For this we need to traverse the tree by choosing the most massive
# progenitor on each level.
#
# While doing this, we can also simultaneously get the number of merger
# events for the most massive halo.
sub find_most_massive {
	my $node = shift;

	# Add the number of mergering events
	push @{$_[1]}, scalar(keys %{$node->[1]});
	
	# Pick the most massive progenitor
	my $max = 0;
	my $maxkey = undef;
	foreach my $key (keys %{$node->[1]}) {
		if ($node->[1]->{$key}->[0][8] > $max) {
			$max = $node->[1]->{$key}->[0][8];
			$maxkey = $key;
		}
	}

	# Return if nothing was found
	return unless (defined $maxkey);

	# Add the index to the table reference given as the first argument.
	push @{$_[0]}, $maxkey;

	# Then recursively continue searching from the most massive progenitor
	&find_most_massive($node->[1]->{$maxkey}, @_);
}

# So, first add the root halo.
my @most_massive = ($halo);
my @mergers;

# Then go though the tree and get the most massive halo at each redshift,
# and add it's mergers
&find_most_massive($tree{$halo}->root_node(), \@most_massive, \@mergers);

# Then get the subhalo data for these most massive halos of their level.
my @subhalos;
for (my $i=0; $i < @most_massive; $i++) {
	push @subhalos, scalar(@{$shdata[$i]->{$most_massive[$i]}});
}


# Now we have all the data to draw a vector graphics image of the merger tree.
# What vector graphics language to use is somewhat arbitrary, but for now, 
# MetaPost is used. Other choices can be easily accommodated.

# Output the commands to start a new merger tree picture
&start_halo($halo);

# Draw the tree.
&draw_tree_branches($tree{$halo}->root_node, 0);

# Then go through level by level, and output the virial radius circles and the
# redshift values.
#&draw_parameter_circles(\@z0_params, \@level_view, $maximum_width);
&draw_parameter_circles($z0_param, \@level_view, $maximum_width);
&draw_labels(\@level_view, \@subhalos, \@mergers, $maximum_width, $p_index);

# Output any needed commands to end this picture
&end_halo($halo);;

}

# Write any possible closing commands
&finalize;

# All halos drawn. Exit successfully.
exit;

#
# The actual drawing is done inside these subroutines. Change these to suit
# other vector graphics languages than MetaPost.
###########################################################################
sub output_preamble {
	print 
	"prologues := 2;\n",
	"u=1cm;\n";
}

# A function to recursively traverse the tree and output the connecting lines
sub draw_tree_branches {
	my $node = shift;
	my $level = shift;

	# Draw a line to every progenitor. 
	foreach my $key (keys %{$node->[1]}) {
		&metapost_line($node->[3], $level, 
			$node->[1]{$key}[3], $level+1);
		&draw_tree_branches($node->[1]{$key}, $level+1);
	}
}

sub draw_parameter_circles {
	my $z0_param = shift;
	my @level_view = @{shift;};
	my $maximum_width = shift;
	for (my $i=0; $i < @level_view; $i++) {
		foreach my $halo (@{$level_view[$i]}) {
			&metapost_circle($halo->[3], $i, 
				$halo->[0][$p_index]/$z0_param);
		}
	}
}


# Draw the labels, using the given parameter
sub draw_labels {
	my $level_view = shift;
	my $subhalos = shift;
	my $mergers = shift;
	my $maximum_width = shift;
	my $p_index = shift;
	for (my $i=0; $i < @$level_view; $i++) {
		# Draw redshift labels
		&metapost_label(-$maximum_width-1.5, $i, 
			sprintf("%.3f", $z_vals[$i]));
	}
	for (my $i=0; $i < @$subhalos; $i++) {
		# Draw subhalo amount labels
		&metapost_label(-$maximum_width-5, $i, 
			sprintf("%d", $subhalos->[$i]));
	}
	for (my $i=0; $i < @$mergers; $i++) {
		# Draw merger event amount labels
		&metapost_label(-$maximum_width-8.5, $i, 
			sprintf("%d", $mergers->[$i]));
	}
	# Then draw a column identifier labels on the top
	&metapost_label(-$maximum_width-1.5, scalar(@$level_view+0.3), 
		"Redshift");
	&metapost_label(-$maximum_width-5, scalar(@$level_view+0.3), 
		"Subhalo");
	&metapost_label(-$maximum_width-5, scalar(@$level_view), 
		"amount");
	&metapost_label(-$maximum_width-8.5, scalar(@$level_view+0.3), 
		"Merger");
	&metapost_label(-$maximum_width-8.5, scalar(@$level_view), 
		"events");
}

sub finalize {
	print "bye\n";
}

sub start_halo {
	my $halo = shift;
	print "beginfig($halo);\n",
	"pickup pencircle scaled 0.001u;\n";
}

sub end_halo {
	my $halo = shift;
	print "endfig;\n";
}



#
# METAPOST drawing primitives
#############################

# Get suitable coordinate scaling and translation
sub scale_coords {
	my $width = shift;
	my $height = shift;

	return ($width*$x_scaling, $height*$y_scaling);
}
sub translate_coords {
	my $width = shift;
	my $height = shift;

	return (6.2-$width, $height+0.5);
}

# Convert from tree coordinates to metapost coordinates
sub convert_coords {
	my $width = shift;
	my $height = shift;

	#return (8.2-$width*0.3, $height*0.9+0.5);
	return &translate_coords(&scale_coords($width, $height));
}

# Draw a line from height0, width0 to height1, width1
sub metapost_line {
	# Convert to xy-coordinates for metapost
	my @start = &convert_coords(@_[0,1]);
	my @end = &convert_coords(@_[2,3]);

	printf "draw (%fu,%fu)--(%fu,%fu);\n", 
		$start[0], $start[1], $end[0], $end[1];
}

sub metapost_circle {
	my @center = &convert_coords(@_[0,1]);
	my $radius = $_[2];

	# Apply the default widthwise coordinate scaling to the radius, and
	# then scale upwards a bit. This is done because the default separation
	# between halos is equal to the sum of their diameters, and the tree
	# looks a bit empty if radius is left as is.
	($radius, undef) = &scale_coords($radius, 0);

	printf "draw fullcircle scaled %fu shifted (%fu,%fu);\n",
		$circle_scaling*$radius, $center[0], $center[1];
}

sub metapost_label {
	my $x = shift;
	my $y = shift;
	my $str = shift; 
	my @coords = &convert_coords($x, $y);

	printf "label.rt(\"$str\", (%fu, %fu));\n", $coords[0], $coords[1];
}
