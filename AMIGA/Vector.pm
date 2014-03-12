# Vector.pm. Written in 2006 by Pauli Pihajoki
package AMIGA::Vector;

use warnings;
use strict;
use Carp;

our(@ISA, @EXPORT, @EXPORT_OK, $VERSION);

use vars qw($VERSION);
$VERSION = 0.10;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(); 

# Overload some primitives.
use overload
'='     =>      sub { $_[0]->copy(); },
'cmp'   =>      sub { 	$_[2] ?
               	     	"$_[1]" cmp $_[0]->to_str() :
               		$_[0]->to_str() cmp "$_[1]" },
'neg'	=>	sub { $_[0]->copy()->neg(); },
'-'     =>      sub { 
			my $c = $_[0]->copy(); $_[2] ?
			$c->neg()->add($_[1]) 
			:
			$c->subtr($_[1]) 
		},
'+'     =>      sub { $_[0]->copy()->add($_[1]); },
'*'     =>      sub { $_[0]->copy()->mult($_[1]); },
'""' 	=> 	sub { $_[0]->to_str(); },
;

#
# Instance and class attributes. These can only be accessed via accessors.
#############################################################

# The class specific attributes, with given default values. 
my %ClassAttributes = (
	# Behaviour control
	verbose => 0
);

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
	# Blessed be.
	bless ($self, $class);

	# Argument should be a list of components.
	$self->{data} = [@_];

	return $self;
}

###########
# Methods #
###########

# Makes a copy of the vector
sub copy {
	my $x = shift;
	my $c = ref($x);
	
	my $self = {};
	bless $self, $c;

	$self->{data} = [@{$x->{data}}];
	return $self;
}

# Multiply by real value
sub mult {
	my $self = shift;
	my $a = shift;

	for (my $i=0; $i < @{$self->{data}}; $i++) {
		$self->{data}[$i] *= $a;
	}

	return $self;
}

# Add a vector to self
sub add {
	my $self = shift;
	my $x = shift;

	if (@{$self->{data}} != @{$x->{data}}) {
		croak "add: dimension mismatch";
	}
	
	for (my $i=0; $i < @{$x->{data}}; $i++) {
		$self->{data}[$i] += $x->{data}[$i];
	}

	return $self;
}

# Subtract a vector from self
sub subtr {
	my $self = shift;
	my $x = shift;

	if (@{$self->{data}} != @{$x->{data}}) {
		croak "subtr: dimension mismatch";
	}
	
	for (my $i=0; $i < @{$x->{data}}; $i++) {
		$self->{data}[$i] -= $x->{data}[$i];
	}

	return $self;
}

# Negate self
sub neg {
	my $self = shift;

	for (my $i=0; $i < @{$self->{data}}; $i++) {
		$self->{data}[$i] = -$self->{data}[$i];
	}

	return $self;
}

# Dot product of two vectors
sub dot {
	my $self = shift;
	my $x = shift;

	if (@{$self->{data}} != @{$x->{data}}) {
		croak "dot: dimension mismatch";
	}

	my $retval = 0;
	for (my $i=0; $i < @{$self->{data}}; $i++) {
		$retval += $self->{data}[$i] * $x->{data}[$i];
	}
	return $retval;
}

# Euclidean norm of vector
sub length {
	my $self = shift;

	return sqrt($self->dot($self));
}

# Create a string representation
sub to_str {
	my $self = shift;

	return "[".join(" ", @{$self->{data}})."]";
}

# Return the dimension of the vector
sub dim {
	my $self = shift;

	return scalar(@{$self->{data}});
}

# Return the components as an array
sub to_array {
	my $self = shift;

	return @{$self->{data}};
}

1;

__END__

=head1 NAME

AMIGA::Vector - Object oriented module for vectors and vector operations.

=head1 SYNOPSIS

	use AMIGA::Vector;

	# Create vectors in 3D.
	my $x = AMIGA::Vector->new(1, 2, 3);
	my $y = AMIGA::Vector->new(-1, 0, 0);

	# Addition and substraction operations are overloaded
	# for ease of use.
	my $z = $x + $y; # $z is now [0, 2, 3]
	print "x - y: ", $x-$y, "\n"; # prints [2 2 3]
	print "-x: ", -$x, "\n"; # prints [-1 -2 -3]

	# They can be called explicitly also, and can be
	# combined
	$z->add($y);   # $z now [-1, 2, 3]
	$z->subtr($y); # $z back to [0, 2, 3]
	$z->add($y)->subtr($y); # leaves $z unchanged

	# Dimensions must match
	$z = AMIGA::Vector->new(1, 2, 3, 4);
	print "x + z", $x + $z, "\n"; # Causes an error!

	# Vector operations
	print $x->dot($y), "\n"; # prints -1
	print $x->length(), "\n"; # prints 3.74165738677394
	print sqrt($x->dot($x)), "\n"; # same as above
	print "dimension: ", $x->dim(), "\n"; # prints 3

	# The components can be retrieved as an array if needed
	my @components = $x->to_array();

=head1 INITIALISATION

=over 4

=item new ( COMPONENTS )

This is a constructor for a new AMIGA::Vector object. C<COMPONENTS> should
be an array of the component values of the vector.

=back

=head1 METHODS

=head2 copy ( )

Returns a copy of this vector.

Arguments: -

Return value: Copy of the vector for which this method was invoked.

=head2 mult ( REAL )

Multiplies this vector with the given scalar C<REAL>.
Note: Changes the vector for which the method was invoked.

Arguments: C<REAL> -- a scalar value.

Return value: The resulting vector.

=head2 add ( VECTOR )

Adds the vector C<VECTOR> to this vector, and then returns it.
Note: Changes the vector for which the method was invoked.

Arguments: C<VECTOR> -- an C<AMIGA::Vector> object.

Return value: The resulting vector.

=head2 subtr ( VECTOR )

Subtracts the vector C<VECTOR> from this vector, and then returns it.
Note: Changes the vector for which the method was invoked.

Arguments: C<VECTOR> -- an C<AMIGA::Vector> object.

Return value: The resulting vector.

=head2 neg ( )

Changes sign of every component of this vector, and also returns it.
Note: Changes the vector for which the method was invoked.

Arguments: -

Return value: The resulting vector.

=head2 dot ( VECTOR )

Calculates the dot product between this vector and the vector
C<VECTOR> given as argument.

Arguments: C<VECTOR> -- an C<AMIGA::Vector> object.

Return value: Dot product of the two vectors.

=head2 length ( )

Calculates the length of this vector and returns it.

Arguments: -

Return value: Length of the vector.

=head2 to_str ( )

Renders the vector into a Perl string. This is done by printing the component
values sequentially separated by whitespace inside enclosing square brackets.
Thus a vector with components C<2.1, 3.2, 4.3> renders into C<[2.1 3.2 4.3]>.

Arguments: -

Return value: Vector string representation.

=head2 dim ( )

Returns the dimension of the vector.

Arguments: -

Return value: Dimension of the vector.

=head2 to_array ( )

Returns the components of the vector as a Perl array.

Arguments: -

Return value: Vector components as an array.

=head1 AUTHOR

Written in 2006 for the benefit of Tuorla Observatory Cosmological Simulations
research group by Pauli Pihajoki.

=head1 COPYRIGHT

Copyright (c) 2006 Pauli Pihajoki. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms
as Perl itself.

=cut

