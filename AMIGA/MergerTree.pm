# MergerTree.pm. Written in 2006 by Pauli Pihajoki
package AMIGA::MergerTree;

use warnings;
use strict;
use Carp;
our(@ISA, @EXPORT, @EXPORT_OK, $VERSION);

use vars qw($VERSION);
$VERSION = 0.10;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(); 

#
# Instance and class attributes. These can only be accessed via accessors.
#############################################################

# The instance specific attributes, with given default values. 
my %InstanceAttributes = (
	# State
	_DOUBLE_LINKED => 0,
	_TREE =>  [],
	# Behaviour control
	verbose => 0
);

# Generate accessors for instance attributes.
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

	# We are a blessed hash reference, as is usual.
	my $self = {};
	# Get a copy of the initial parameters.
	%$self = %InstanceAttributes;
	# Blessed be.
	bless ($self, $class);

	# The only argument should be an array reference containing the
	# data for the root node and a hash reference of the progenitors.
	if (@_ != 1) {
		confess "new: arguments: array table reference";
	}
	my $ref = shift;;

	# Get the data
	$self->{_TREE} = $ref;

	return $self;
}

###########
# Methods #
###########

# Return the first node of the tree.
sub root_node {
	my $self = shift;
	return $self->{_TREE};
}

# Make the tree doubly linked. Note: this is irreversible.
sub double_link {
	my $self = shift;
	
	# Don't double link if we already are.
	return if ($self->{_DOUBLE_LINKED});

	# First, add an empty hash reference for the root node
	push @{$self->{_TREE}}, {};

	# Then update all progenitors.
	&create_links($self->root_node());

	# We are now double linked
	$self->{_DOUBLE_LINKED} = 1;
}

# Create and return a by-depth-level view of the tree.
sub level_view {
	my $self = shift;

	# Start with the root node at the first level.
	my @level_view = ([$self->root_node()]);

	# Then go through the tree recursively and update
	&update_level_view($self->root_node(), \@level_view, 0);

	# Return the result;
	return @level_view;
}

# Return the number of branchings from the given node, or from the root
# node, if no arguments.
sub branchings {
	my $self = shift;
	my $node;
	if (@_ > 1) {
		$node = shift;
	}
	else {
		$node = $self->root_node;
	}

	# Go through the tree recursively, incrementing the branching counter,
	# and return the result.
	return &tree_size($node);
}

# Map some function to each node of the tree by doing a recursive traversal
sub traverse_map {
	my $self = shift;
	my $func = shift;
	my @args = @_;

	# Apply the function recursively, return the result
	return &$func($self->root_node(), @args);
}

####################
# Helper functions #
####################

# For double linking
sub create_links {
	my $node = shift;

	# Foreach progenitor, add a reference to this node. Then process
	# each progenitor in the same way
	foreach my $key (keys %{$node->[1]}) {
		push @{$node->[1]{$key}}, $node; 
		&create_links($node->[1]{$key});
	}
}

# For making a by-level-view
sub update_level_view {
	my $node = shift;
	my $view = shift;
	my $level = shift;

	# If there is nothing yet on the next level, and we will add something,
	# then push an new array reference
	if ($#$view <= $level && (keys %{$node->[1]}) != 0) {
		push @$view, [];
	}

	# Add every progenitor node, and also add the index numbers
	foreach my $key (keys %{$node->[1]}) {
		push @{$view->[$level+1]}, $node->[1]{$key};
		&update_level_view($node->[1]{$key}, $view, $level+1);
	}
}

# A function to calculate the number of downstream branches from a given
# node.
sub tree_size {
	my $node = shift;

	# If there are no progenitors, terminate with 1.
	if ((keys %{$node->[1]}) <= 0) {
		return 1;
	}

	# Otherwise, return the sum of progenitor sizes.
	my $result = 0;
	foreach my $key (keys %{$node->[1]}) {
		$result += &tree_size($node->[1]{$key});
	}

	return $result;
}


1;

__END__

=head1 NAME

AMIGA::MergerTree - An object oriented tree structure representation of
halo data.

=head1 SYNOPSIS

	use AMIGA::MergerTree;

	my $tree = AMIGA::MergerTree->new(\%my_tree_structure);

=head1 INITIALISATION 

=over 4

=item new ( TREE )

This is the constructor for a new AMIGA::MergerTree object. C<TREE> is 
a reference to a hash table containing the recursive tree structure.

This structure should consist of a hash table with halo indexes pointing to a
array references containing:

=over 4

=item *

Array reference to the _halos data for the halo.

=item *

A hash table with progenitor halo indexes pointing to array references
containing:

=over 4

=item *

Array reference to the _halos data for the halo.

=item *

A hash table with the progenitor indexes of this halo pointing to identical
array references as above, with the recursion continuing through the tree.

=back

=back

The tree structure of the MergerTree is so an unordered tree, that can be
only walked depthwards, i.e. from root node towards leaf nodes.

=back

=head1 METHODS

=head2 root_node ( )

Returns the first node in the raw tree structure, as described in 
L</INITIALISATION>.

Arguments: -

Return value: The root node of the tree structure.

=head2 double_link ( )

Makes the whole tree doubly linked, adding to each node array a hash reference
of the parent node. This allows walking the tree in both directions.

Note: This method alters the object that it was invoked for. The process is
irreversible.

Arguments: -

Return value: -

=head2 level_view ( )

Returns a 2D array of the contents of the tree structure, with
the first dimension being the depth level of the tree, and the second the
nodes at that depth. Thus the number of columns at each row equals the number
of nodes at that depth.

Arguments: -

Return value: A 2D array of tree contents.

=head2 branchings ( [NODE] )

Returns the total number of branchings that occur in the tree after the node
given in C<NODE> or after the root node, if called with no arguments.

Arguments: C<NODE> -- reference to an array containing node information. See
L</INITIALISATION>.

Return value: Total number of branchings depthwards from the given node
or from the root node.

=head1 AUTHOR

Written in 2006 for the benefit of Tuorla Observatory Cosmological Simulations
research group by Pauli Pihajoki.

=head1 COPYRIGHT

Copyright (c) 2006 Pauli Pihajoki. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms
as Perl itself.

=cut

