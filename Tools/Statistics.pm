# Statistics.pm. Written in 2006 by Pauli Pihajoki
# Functionality for basic statistic operations, and an interface for
# arbitrary row and column operations.
package Tools::Statistics;

use warnings;
use strict;
use Carp;

our(@ISA, @EXPORT, @EXPORT_OK, $VERSION);

use vars qw($VERSION);
$VERSION = 0.01;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(mean median variance std_deviation kahan_sum); 
@EXPORT_OK = qw();

# Global constants
my %Constants = (
);

# Generate accessors for global constants
for my $datum (keys %Constants) { 
	no strict "refs";    
	*$datum = sub {
		use strict "refs";    
		my ($newvalue) = @_;
		$Constants{$datum} = $newvalue if @_ > 0;
		return $Constants{$datum};
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

	# Argument should be a 2D array. Go through it and make sure it's
	# valid.
	croak "new: arguments: 2D array" if (@_ == 0 || @{$_[0]} == 0 
		|| ref($_[0]->[0]));
	my $cols;
	$cols = scalar(@{$_[0]});
	for my $row (@_) {
		croak "new: argument not a valid 2D array" 
			unless (ref($row) eq "ARRAY");
		croak "new: found only ",scalar(@$row)," columns, expected ",
		      $cols, " columns" if (@$row != $cols);
	}

	# Looks ok. Set the data, and create column view.
	$self->{_DATA} = [@_];
	$self->{_COLS} = [];
	for (my $i=0; $i < @_; $i++) {
		for (my $j=0; $j < $cols; $j++) {
			push @{$self->{_COLS}->[$j]}, $self->{_DATA}->[$i][$j];
		}
	}


	return $self;
}

# Accessors for amount of rows and columns
sub rows {
	my $self = shift;

	return scalar(@{$self->{_DATA}});
}

sub cols {
	my $self = shift;

	return scalar(@{$self->{_COLS}});
}

# Perform some arbitrary operation for each given column, or for all columns
# if not specified.
# Return array of results.
sub map_to_columns {
	my $self = shift;
	croak "map_to_columns: arguments: operation [, columns]" if (@_ < 1);

	my $op = shift;

	my @cols;
	if (@_) {
		@cols = @_;
	}
	else {
		@cols = (0..($self->cols()-1));
	}

	# For each column, perform the operation and return result
	my @result;
	for (@cols) {
		croak "map_to_columns: column index $_ too large (max. ",
		 scalar(@{$self->{_COLS}})-1, ")" if ($_ >= @{$self->{_COLS}});
		push @result, &$op(@{$self->{_COLS}->[$_]});
	}

	return @result;
}

# Perform some arbitrary operation for given rows or each row if not
# specified.
# Return array of results.
sub map_to_rows {
	my $self = shift;
	croak "map_to_rows: arguments: operation [, rows]" if (@_ < 1);

	my $op = shift;

	my @rows;
	if (@_) {
		@rows= @_;
	}
	else {
		@rows = (0..($self->rows()-1));
	}

	# Operate on the rows, and collect the results
	my @result;
	for (@{$self->{_DATA}}[@rows]) {
		push @result, &$op(@$_);
	}

	return @result;
}

#
# Some predefined statistical operations
########################################

# Mean for given array
sub mean {
	my @data = @_;

	# Sort ascending, sum and then divide.
	@data = sort {$a <=> $b} @data;
	#my $sum = 0;
	#for (@data) {
	#	$sum += $_;
	#}
	my $sum = &kahan_sum(@data);

	# Force floating point evaluation.
	return $sum/(1.0*@data);
}

# Median
sub median {
	my @data = @_;

	# Sort ascending and return result
	@data = sort {$a <=> $b} @data;

	if (@data % 2) {
		return $data[@data/2];
	}

	return 0.5*($data[@data/2 - 1] + $data[@data/2]);
}

# Variance 
sub variance {
	my @data = @_;

	# Get mean and calculate variance.
	my $mean = &mean(@data);

	# Get the square of differences between value and mean
	for (@data) {
		$_ = ($_ - $mean)**2;
	}
	# Sum them
	my $sum = &kahan_sum(@data);

	# Return result
	return $sum/@data;
}

# Standard deviation
sub std_deviation {
	return sqrt(&variance(@_));
}

# Kahan compensated summing algorithm
sub kahan_sum {
	croak "kahan_sum: no input" unless (@_);
	my $c = 0.0;
	my $sum = $_[0];
	my ($y, $t);
	for (my $i=1; $i<@_; $i++) {
		#printf "c: %16.16e\n", $c;
		$y = $_[$i] - $c;
		#printf "sum: %16.16e y: %16.16e\n", $sum, $y;
		$t = $sum + $y;
		#printf "t: %16.16e\n", $t;
		# Order of calculation is very important here
		$c = $t - $sum;
		#printf "c: %16.16e\n", $c;
		$c -= $y;
		#printf "c: %16.16e\n", $c;
		$sum = $t;
	}

	return $sum;
}






1;

__END__

=head1 NAME

Tools::Statistics - Simple object oriented statistics module for handling 2D
data formats.

=head1 SYNOPSIS

	use Tools::Statistics;

	my $stat = Tools::Statistics->new(@some_2d_data);

	# Find out number of rows and columns
	my $num_rows = $stat->rows();
	my $num_cols = $stat->cols();

	# Access basic statistics functions on a per-column basis
	my @column_means = $stat->map_to_columns(\&mean);
	my @column_medians = $stat->map_to_columns(\&median);
	my @col_variances = $stat->map_to_columns(\&variance);
	my @col_std_devs = $stat->map_to_columns(\&std_deviation);

	# Or just for some columns
	my @col_0_1_means = $stat->map_to_columns(\&mean, 0, 1);

	# Same things can be done for rows
	my @row_means = $stat->map_to_rows(\&mean);
	my @row_0_1_means = $stat->map_to_rows(\&mean, 0, 1);

	# Apply arbitrary operations on a per-column or per-row basis
	my $sum_of_args = sub { 
		my $ret=0; 
		for (@_) {
			$ret += $_;
		} 
		return $ret; 
	};
	my @column_sums = $stat->map_to_column($sum_of_args);
	my @row_sums = $stat->map_to_rows($sum_of_args);

	# Apply some arbitrary operation between different columns on 
	# a row by row basis and collect the results
	my @division_of_col_3_with_col_4 = 
			$stat->map_to_rows( { $_[3] / $_[4]; } );

=head1 INITIALISATION 

=over 4

=item new ( DATA )

This is a constructor for a new Tools::Statistics object. C<DATA> is
a 2D array (with uniform amount of columns on each row).

The constructor does some primitive validation on the data, and creates
a column-based view on the data.

=back

=head1 METHODS

=head2 rows ( )

Returns the number of rows in the data.

Arguments: -

Return value: Number of rows in the data.

=head2 cols ( )

Returns the number of columns in the data.

Arguments: -

Return value: Number of columns in the data.

=head2 map_to_columns ( OP [, COLUMN [, COLUMN [,...]]] )

Maps some operation C<OP> to every element in given columns, or all columns
if none are specified. The results are gathered in a list and returned.

Arguments:

=over 4

=item *

C<OP> -- A subroutine reference. The subroutine should expect to have an 
array of all the elements of a single column passed as an argument. These
will be in the same order as in the original data, with topmost in 
C<$_[0]>. The return value of the subroutine will be stored in the result
array.

=item *

C<COLUMN> -- A number specifying the column index. Indexing starts from 0.

=back

Return value: A list of the results of all the evaluations of C<OP>. These
results are in the same order as the given column indexes or in the same
order (left to right) as in the original data, if no columns were specified.

=head2 map_to_rows ( OP [, ROW [, ROW [,...]]] )

Maps some operation C<OP> to every element in given rows, or all rows 
if none are specified. The results are gathered in a list and returned.

Arguments:

=over 4

=item *

C<OP> -- A subroutine reference. The subroutine should expect to have an 
array of all the elements of a single row passed as an argument. These
will be in the same order as in the original data, with leftmost in 
C<$_[0]>. The return value of the subroutine will be stored in the result
array.

=item *

C<COLUMN> -- A number specifying the column index. Indexing starts from 0.

=back

Return value: A list of the results of all the evaluations of C<OP>. These
results are in the same order as the given row indexes, or in the same
order as in the original data, if no rows were specified.

=head1 UTILITY FUNCTIONS

These take input data as an array of values (C<VALUES>) and return the result
of the corresponding statistical operation on that data. All summing is
done by sorting to ascending order and using the compensating Kahan summation
algorithm to try to minimise summation errors.

These can be mapped to columns or rows as described in L</SYNOPSIS>.

=head2 mean ( VALUES )

Calculates the mean of C<VALUES>.

=head2 median ( VALUES ) 

Calculates the median of C<VALUES>.

=head2 variance ( VALUES )

Calculates the variance of C<VALUES>. Note: The denominator is taken as C<N>,
and not C<N-1> so this is not the unbiased sample variance.

=head2 std_deviation ( VALUES )

Calculates the standard deviation of C<VALUES>. Note: The denominator is taken
as C<N>, and not C<N-1> so this is not the unbiased sample standard deviation.

=head1 AUTHOR

Written in 2006 for the benefit of Tuorla Observatory Cosmological Simulations
research group by Pauli Pihajoki.

=head1 COPYRIGHT

Copyright (c) 2006 Pauli Pihajoki. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms
as Perl itself.

=cut

