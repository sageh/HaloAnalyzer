# Random.pm. Written in 2006 by Pauli Pihajoki
package Tools::Random;

use warnings;
use strict;
use Carp;
#use Math::TrulyRandom;

our(@ISA, @EXPORT, @EXPORT_OK, $VERSION);

use vars qw($VERSION);
$VERSION = 0.01;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(); 

#
# Instance and class attributes. These can only be accessed via accessors.
#############################################################

# The instance specific attributes, with given default values. 
my %InstanceAttributes = (
	# Behaviour control
	verbose => 0
);

# Class specific attributes
my %ClassAttributes = (
	# The random seed common for every instance
	_SEED => undef,
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

# Generate accessors for class attributes.
for my $datum (keys %ClassAttributes) { 
	no strict "refs";    
	*$datum = sub {
		use strict "refs";    
		my ($self, $newvalue) = @_;
		$ClassAttributes{$datum} = $newvalue if @_ > 1;
		return $ClassAttributes{$datum};
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

	# If the RNG hasn't been seeded yet, then seed it.
	unless (defined $self->_SEED) {
		print "Generating random seed...\n" if $self->verbose;
		$self->_SEED(time ^ $$ ^ unpack "%L*", `ps axww | gzip`);
		srand $self->_SEED;
	}

	return $self;
}

###########
# Methods #
###########

# Depending on the arguments, return a random real number between 0 and the
# argument, or between the arguments.
#
# Arguments: one or two real numbers
sub rand {
	my $self = shift;

	if (@_ >= 2) {
		return $_[0] + rand ($_[1] - $_[0]);
	}
	elsif (@_ == 1) {
		return rand $_[0];

	}

	# If no arguments, work like perl rand
	return rand;
}

# Returns an d-dimensional binomial random field of N points, with each
# component ranging between the given interval. Returns an two dimensional
# array indexed by points and vector components.
#
# Arguments: hash table:
#   
#   dimension => number of dimensions (d) (default: 1)
#   points => number of points (N) 	  (default: 1)
#   low_bound => [lb_1, lb_2, ..., lb_d]  (default: [0, 0, ..., 0])
#   high_bound => [hb_1, hb_2, ..., hb_d] (default: [1, 1, ..., 1])
sub binomial_field {
	my $self = shift;

	if (@_ % 2) {
		croak "binomial_field: arguments: hash table";
	}

	my %args = @_;

	unless (exists $args{points}) {
		$args{points} = 1;
	}

	unless (exists $args{dimension}) {
		$args{dimension} = 1;
	}

	if (exists $args{low_bound}) {
		if (@{$args{low_bound}} != $args{dimension}) {
			croak "binomial_field: low bound has ".
			scalar(@{$args{low_bound}})." components, but there ".
			"are $args{dimension} dimensions";
		}
	}
	else {
		# Default low bound is all zeroes.
		$args{low_bound} = [];
		for (my $i=0; $i < $args{dimension}; $i++) {
			push @{$args{low_bound}}, 0;
		}
	}

	if (exists $args{high_bound}) {
		if (@{$args{high_bound}} != $args{dimension}) {
			croak "binomial_field: high bound has ".
			scalar(@{$args{high_bound}})." components, but there ".
			"are $args{dimension} dimensions";
		}
	}
	else {
		# Default high bound is all ones.
		$args{high_bound} = [];
		for (my $i=0; $i < $args{dimension}; $i++) {
			push @{$args{high_bound}}, 1;
		}
	}

	# Build the array.
	my @result;

	print "Generating binomial field...\n" if $self->verbose;
	for (my $i=0; $i < $args{points}; $i++) {
		push @result, [];
		for (my $j=0; $j < $args{dimension}; $j++) {
			push @{$result[$#result]}, $self->rand(
				$args{low_bound}->[$j],
				$args{high_bound}->[$j]);
		}
	}

	return @result;
}

1;

__END__

=head1 NAME

Tools::Random - Simple object oriented module for random number functionality.

=head1 SYNOPSIS

	use Tools::Random;

	# Create a new randomizer. This will also seed the random
	# number generator (RNG) if it hasn't been seeded yet.
	my $rnd = Tools::Random->new(verbose => 1);

	# A random number from range [$x, $y]
	my $xynum = $rnd->rand($x, $y);

	# A random number from range [0, $x]
	my $xnum = $rnd->rand($x);

	# Works like a call to perl rand()
	my $num = $rnd->rand();

	# Create a binomial random field of 1000 points in 3 dimensions,
	# (pseudo)uniformly distributed between given limits. Return
	# it in a 2D array.
	my @field = $rnd->binomial_field(
				points => 1000,
				dimensions => 3,
				low_bound => [-5, -3, -4],
				high_bound => [3, 4, 5]);

=head1 INITIALISATION

=over 4

=item new ( OPTIONS )

This is a constructor for a new Tools::Random object. C<OPTIONS> should be
a hash table with the following possible keys:

=over 4

=item verbose

Controls the verbosity level of the module. Currently only two values are
supported. Value 1 causes the module to output some information of the ongoing
processes to standard output, while value 0 suppresses most output.

Default: 0 

=back

If the Perl RNG has not been seeded yet, the constructor provides a seed
by calling C<srand(time ^ $$ ^ unpack "%L*", `ps axww | gzip`)> which should be
sufficient for low criticality random data generation. This seed is then used
by all instances of the C<Tools::Random> class.

=back

=head1 METHODS

=head2 rand ( [LOW [, HIGH] )

Returns a random floating point number from a I<half-closed interval>, i.e.
for low bound C<x> and high bound C<y> a number from C<[x,y)> is returned.

When both parameters are provided, returns a random number between
C<LOW> and C<HIGH>. If only one parameter is given, a random number between
0 and the parameter is returned. If no parameters are given, a number between
0 and 1 is returned.

Arguments: 

=over 4

=item *

C<LOW> -- random number low bound, inclusive.

=item *

C<HIGH> -- random number high bound, exclusive.

=back

Return value: A random floating point value between C<LOW> and C<HIGH>, 0 and
C<HIGH> or 0 and 1.

=head2 binomial_field ( OPTIONS )

Creates a series of random floating point vectors from an N-dimensional space
between given limits. These are returned as a 2D Perl array.

This method uses the standard Perl RNG to generate the data, and depending
on the local implementation and the criticality of the application could 
produce unsatisfactory data. It is recommended that a separate method using
for example the Mersenne Twister algorithm be used for applications where
high quality pseudorandom numbers are required.

Arguments:

=over 4

=item *

C<OPTIONS> -- a hash table with following possible keys:

=over 4

=item dimension

The dimension N of each random vector.

Default: 1

=item low_bound

The low limiting bound of the random vectors, inclusive. 
Given as an array reference with N components.

Default: [0]

=item high_bound

The high limiting bound of the random vectors, exclusive. 
Given as an array reference with N components.

Default: [1]

=item points

The number of vectors to create.

Default: 1

=back

=back

=head1 AUTHOR

Written in 2006 for the benefit of Tuorla Observatory Cosmological Simulations
research group by Pauli Pihajoki.

=head1 COPYRIGHT

Copyright (c) 2006 Pauli Pihajoki. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms
as Perl itself.

=cut

