# Correlation.pm. Written in 2006 by Pauli Pihajoki
package AMIGA::Correlation;

use warnings;
use strict;
use Carp;
use Tools::Random;
use Tools::Vector;

our(@ISA, @EXPORT, @EXPORT_OK, $VERSION);

use vars qw($VERSION);
$VERSION = 0.10;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(); 
@EXPORT_OK = qw(calculate_correlation);

# Global constants
my %Constants = (
	# How many more particles should the automatically generated comparison
	# random binomial field have than the data sea.
	# XXX: 10 is apparently a decent default value, but DR correlations
	# could benefit from more.
	#random_field_multiplier => 10,
	random_field_multiplier => 20,
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

###############
# Module code #
###############

# Calculate a Landy-Szalay (1993) correlation for given data. Returns an array
# with the given bins and the correlation function value for that bin.
#
# Arguments: hash with:
# low_bound => data low bound (Vector)
# high_bound => data high bound (Vector)
# data => data [N] (Vector)
# bins => (parallel) scale bins[x][2] 
# field => (binomial field (arrayref))

sub calculate_correlation {
	if (@_ % 2) {
		croak "calculate_correlation: argument isn't reducible ",
		"to a hash";
	}
	my %args = @_;

	# Get the arguments
	my $low_bound = $args{low_bound};
	my $high_bound = $args{high_bound};
	my $data = $args{data};
	my $bins = $args{bins};

	my $trbins = undef;
	if (exists $args{trbins}) {
		if (defined $args{trbins}) {
			$trbins = $args{trbins};
		}
		else {
			carp "calculate_correlation: transverse bins ",
			"specified, but value is undefined";
		}
	}
	my $field = undef;
	if (exists $args{field}) {
		if (defined $args{field}) {
			$field = $args{field};
		}
		else {
			carp "calculate_correlation: field specified, ",
			"but value is undefined";
		}
	}
	my $ext_prg = undef;
	if (exists $args{external_prg}) {
		if (defined $args{external_prg}) {
			$ext_prg = $args{external_prg};
		}
		else {
			carp "calculate_correlation: external correlator ",
		        "specified, but value is undefined";
		}
	}

	# Store the number of data points separately for convenience
	my $N = @$data;
	
	if ($N <= 0) {
		carp "calculate_correlation: no data";
		return undef;
	}
	if (@$bins <= 0) {
		carp "calculate_correlation: no bins";
		return undef;
	}
	if (defined $trbins && @$trbins <= 0) {
		carp "calculate_correlation: no bins";
		return undef;
	}
	if ($low_bound->dim() != $high_bound->dim()
		|| $low_bound->dim() != $data->[0]->dim()) {
		carp "calculate_correlation: dimension mismatch";
		return undef;
	}

	# We calculate the value for each bin. If the bin is defined by:
	# x \in [b_low, b_high), then the value we return is ksi(b_low).
	# This means that bins can for example be given in the form [r, dr).
	# Not that the bin includes lower bound, but excludes higher bound.

	# Use this many points for the comparison random data
	my $N_r = $N * &random_field_multiplier(); 

	# Generate the binomial random field for R-pairs, if not given one
	my @r;
	if (defined $field) {
		@r = @$field;
	}
	else {
		my $random = Tools::Random->new(verbose=>1);
		@r = $random->binomial_field(
			dimension => $low_bound->dim(),
			points => $N_r,
			low_bound => $low_bound->{data},
			high_bound => $high_bound->{data});
	}

	# Vectorize the field for easier use, as the data that's passed in
	# is already expressed via the Vector class
	for (my $i=0; $i < @r; $i++) {
		$r[$i] = Tools::Vector->new(@{$r[$i]});
	}

	# If we want to use the external implementation for calculating
	# the pairs and the correlation function, then diverge to it here
	if (defined $ext_prg) {
		my @ext; 
		if (defined $trbins) {
			@ext = &external_correlator($ext_prg, 
			$bins, $trbins, $low_bound, $high_bound, $data, \@r);
		}
		else {
			@ext = &external_correlator($ext_prg, 
			$bins, $low_bound, $high_bound, $data, \@r);
		}
		return @ext;
	}

	# Otherwise, proceed with using perl for the calculation.
	if (defined $trbins) {
		return &perl_component_correlator(
			$bins, $trbins, $low_bound, $high_bound, 
			$data, \@r);
	}
	else {
		return &perl_correlator(
			$bins, $low_bound, $high_bound, $data, \@r);
	}
}

sub external_correlator {
	croak "external_correlator: arguments: binary, (parallel) bins, ",
	"[transverse bins], low bound, high bound, data, comparison field" 
		if (@_ != 6 && @_ != 7);

	my ($prg, $bins, $trbins, $lb, $hb, $data, $cmp);
	if (@_ == 6) {
		($prg, $bins, $lb, $hb, $data, $cmp) = @_;
		$trbins = undef;
	}
	else {
		($prg, $bins, $trbins, $lb, $hb, $data, $cmp) = @_;
	}

	# Create suitable temporary files and run the program.
	#$File::Temp::DEBUG = 1;
	my $datafile = new File::Temp(
		TEMPLATE => 'correlator_input_data_XXXXX');
	my $cmpfile = new File::Temp(
		TEMPLATE => 'correlator_input_cmp_XXXXX');
	my $outfile = new File::Temp(
		TEMPLATE => 'correlator_output_XXXXX');

	# Write the number of bins, and the bins.
	# If we were given transverse bins, then write both sets of bins
	print $datafile scalar(@$bins), "\n";
	foreach my $bin (@$bins) {
		print $datafile "$bin->[0] $bin->[1]\n";
	}
	
	if (defined $trbins) {
		print $datafile scalar(@$trbins), "\n";
		foreach my $bin (@$trbins) {
			print $datafile "$bin->[0] $bin->[1]\n";
		}
	}

	# Number of data points, bounds, and data points themselves
	print $datafile scalar(@$data), "\n";
	print $datafile join(" ", @{$lb->{data}}), "\n";
	print $datafile join(" ", @{$hb->{data}}), "\n";
	foreach my $dp (@$data) {
		print $datafile join(" ", @{$dp->{data}}), "\n";
	}

	# Print comparison field

	# Number of data points, bounds, and data points themselves
	print $cmpfile scalar(@$cmp), "\n";
	print $cmpfile join(" ", @{$lb->{data}}), "\n";
	print $cmpfile join(" ", @{$hb->{data}}), "\n";
	foreach my $dp (@$cmp) {
		print $cmpfile join(" ", @{$dp->{data}}), "\n";
	}

	# DEBUG:
	# Try to read what we just wrote.
	#print "$datafile: \n";
	#seek($datafile, 0, 0);
	#while (<$datafile>) {print;}

	# Run the correlator binary, and store the data in the temporary file
	my $runcmd = $prg." $outfile $datafile $cmpfile";
	if (system($runcmd)) {
		carp "Running external correlator $prg failed!";
		return undef;
	}

	# Now collect the output 
	my @result;
	my @tmp;
	seek($outfile, 0, 0);
	while (<$outfile>) {
		chomp;
		next if (/^#/ || /^\s*$/);
		@tmp = split;
		push @result, [@tmp];
	}

	return @result;
}

sub perl_correlator {
	croak "perl_correlator: arguments: bins, low bound, high bound, ",
	"data, comparison field" if (@_ != 5);

	my ($bins, $lb, $hb, $data, $r) = @_;

	my @dd;
	my @dr;
	my @rr;

	# Fill the pair containers with zeroes
	for (@$bins) {
		push @dd, 0;
		push @dr, 0;
		push @rr, 0;
	}

	my $N = scalar(@$data);
	my $N_r = scalar(@$r);

	# Go through each particle in the data set and calculate DD and DR
	# pairs. 
	my $tmpv;
	my $tmpl;
	print "DD and DR pairs: \n";
	OUTER: for (my $i=0; $i < $N; $i++) {
		print "$i ";
		print "\n" if ($i != 0 && ($i % 10 == 0));
		# Calculate DR pairs
		for (my $j=0; $j < @$r; $j++) {
			$tmpv = $data->[$i] - $r->[$j];
			$tmpl = $tmpv->length();
			for (my $k=0; $k < @$bins; $k++) {
				# Check if the pair belongs into this bin
				if ($bins->[$k][0] <= $tmpl
					&& $tmpl <= $bins->[$k][1]) {
					# Add one DR pair for this bin
					$dr[$k]++;
				}
			}
		}

		# Calculate DD pairs
		for (my $j=$i+1; $j < $N; $j++) { 	
			$tmpv = $data->[$i] - $data->[$j];
			$tmpl = $tmpv->length();
			for (my $k=0; $k < @$bins; $k++) {
				# Check if the pair belongs into this
				# bin
				if ($bins->[$k][0] <= $tmpl
					&& $tmpl <= $bins->[$k][1]) {
					# Add one DD pair for this bin
					$dd[$k]++;
				}
			}
		}
	}
	print "\n";

	# Calculate RR pairs
	print "RR pairs: \n";
	OUTER: for (my $i=0; $i < @$r; $i++) {
		# DEBUG:
		print "$i ";
		print "\n" if ($i != 0 && ($i % 10 == 0));
		for (my $j=$i+1; $j < @$r; $j++) {
			$tmpv = $r->[$i] - $r->[$j];
			$tmpl = $tmpv->length();
			for (my $k=0; $k < @$bins; $k++) {
				# Check if the pair belongs into this
				# bin
				if ($bins->[$k][0] <= $tmpl
					&& $tmpl <= $bins->[$k][1]) {
					# Add one RR pair for this bin
					$rr[$k]++;
				}
			}
		}
	}
	print "\n";

	# When self-pairs aren't counted, there are
	# n(n-1) non-unique pairings, and thus for k unique pairs, there
	# are 2*k non-unique pairs, so we just need to multiply by 2.
	
	print "Uncorrected pairs:\n";
	# Print pairs
	for (my $i=0; $i < @$bins; $i++) {
		printf "Bin %d: DD: %d DR: %d: RR: %d\n",
		$i, $dd[$i], $dr[$i], $rr[$i];

		# For non-unique pairs excluding self-pairs
		$dd[$i] *= 2;
		$rr[$i] *= 2;
	}

	# DEBUG:
	# Print pairs
	print "Corrected pairs:\n";
	for (my $i=0; $i < @$bins; $i++) {
		printf "Bin %d: DD: %d DR: %d: RR: %d\n",
		$i, $dd[$i], $dr[$i], $rr[$i];
	}


	# Calculate and store the value of the correlation function, and it's
	# Poisson error (due to Statistics of Galaxy Distribution, p. 91, by 
	# Martinez & Saar) for every bin
	my @result;
	# DEBUG: also print some data
	my ($ksi, $error);
	for (my $i=0; $i < @$bins; $i++) {
		if ($rr[$i] != 0 && $dd[$i] != 0) {
			$ksi = 1 + ($N_r/$N)**2*$dd[$i]/$rr[$i] - 
				2*($N_r/$N)*$dr[$i]/$rr[$i];
			$error = ($N/$N_r)**2 * (1 + $ksi)/sqrt($rr[$i]);
		}
		else {
			carp "RR($bins->[$i][0]) or DD($bins->[$i][0]) was 0";
			$ksi = $error = 0;
			
		}
		push @result, [@{$bins->[$i]}, $ksi, $error];
	}

	return @result;
}

sub perl_component_correlator {
	croak "perl_component_correlator: arguments: bins, trbins, low bound, ",
	"high bound, data, comparison field" if (@_ != 6);

	my ($bins, $trbins, $lb, $hb, $data, $r) = @_;

	my @dd;
	my @dr;
	my @rr;

	# Fill the pair containers with zeroes
	for (@$bins) {
		push @dd, [(0) x scalar(@$trbins)];
		push @dr, [(0) x scalar(@$trbins)];
		push @rr, [(0) x scalar(@$trbins)];
	}

	my $N = scalar(@$data);
	my $N_r = scalar(@$r);

	# Go through each particle in the data set and calculate DD and DR
	# pairs. 
	my ($sep, $los);
	my ($para, $trans);
	print "DD and DR pairs: \n";
	OUTER: for (my $i=0; $i < @$data; $i++) {
		print "$i ";
		print "\n" if ($i != 0 && ($i % 10 == 0));
		# Calculate DR pairs
		for (my $j=0; $j < @$r; $j++) {
			$sep = $data->[$i] - $r->[$j];
			$los = 0.5*($data->[$i] + $r->[$j]);
			$para = $sep->dot($los)/$sep->length();
			$trans = $sep->dot($sep) - $para*$para;
			# Check if the pair belongs into any bins
			DR: for (my $k=0; $k < @$bins; $k++) {
				for (my $l=0; $l < @$trbins; $l++) {
					if ($bins->[$k][0] <= $para &&
						$para <= $bins->[$k][1] &&
						$trbins->[$l][0] <= $trans
						&& $trans <= $trbins->[$l][1]) {
						$dr[$k][$l]++;
						last DR;
					}
				}
			}
		}

		# Calculate DD pairs
		for (my $j=$i+1; $j < @$data; $j++) { 	
			$sep = $data->[$i] - $data->[$j];
			$los = 0.5*($data->[$i] + $data->[$j]);
			$para = $sep->dot($los)/$sep->length();
			$trans = $sep->dot($sep) - $para*$para;
			# Check if the pair belongs into any bins
			DD: for (my $k=0; $k < @$bins; $k++) {
				for (my $l=0; $l < @$trbins; $l++) {
					if ($bins->[$k][0] <= $para &&
						$para <= $bins->[$k][1] &&
						$trbins->[$l][0] <= $trans
						&& $trans <= $trbins->[$l][1]) {
						$dd[$k][$l]++;
						last DD;
					}
				}
			}
		}
	}
	print "\n";

	# Calculate RR pairs
	print "RR pairs: \n";
	OUTER: for (my $i=0; $i < @$r; $i++) {
		# DEBUG:
		print "$i ";
		print "\n" if ($i != 0 && ($i % 10 == 0));
		for (my $j=$i+1; $j < @$r; $j++) {
			$sep = $r->[$i] - $r->[$j];
			$los = 0.5*($r->[$i] + $r->[$j]);
			$para = $sep->dot($los)/$sep->length();
			$trans = $sep->dot($sep) - $para*$para;
			# Check if the pair belongs into any bins
			RR: for (my $k=0; $k < @$bins; $k++) {
				for (my $l=0; $l < @$trbins; $l++) {
					if ($bins->[$k][0] <= $para &&
						$para <= $bins->[$k][1] &&
						$trbins->[$l][0] <= $trans
						&& $trans <= $trbins->[$l][1]) {
						$rr[$k][$l]++;
						last RR;
					}
				}
			}
		}
	}
	print "\n";

	# When self-pairs aren't counted, there are
	# n(n-1) non-unique pairings, and thus for k unique pairs, there
	# are 2*k non-unique pairs, so we just need to multiply by 2.
	
	print "Uncorrected pairs:\n";
	# Print pairs
	for (my $i=0; $i < @$bins; $i++) {
		printf "Bin %d: DD: %d DR: %d: RR: %d\n",
		$i, $dd[$i], $dr[$i], $rr[$i];

		# For non-unique pairs excluding self-pairs
		$dd[$i] *= 2;
		$rr[$i] *= 2;
	}

	# DEBUG:
	# Print pairs
	print "Corrected pairs:\n";
	for (my $i=0; $i < @$bins; $i++) {
		printf "Bin %d: DD: %d DR: %d: RR: %d\n",
		$i, $dd[$i], $dr[$i], $rr[$i];
	}


	# Calculate and store the value of the correlation function, and it's
	# Poisson error (due to Statistics of Galaxy Distribution, p. 91, by 
	# Martinez & Saar) for every bin
	my @result;
	my ($ksi, $error);
	for (my $i=0; $i < @$bins; $i++) {
		for (my $j=0; $j < @$trbins; $j++) {
			if ($rr[$i] != 0 && $dd[$i] != 0) {
				$ksi = 1 + ($N_r/$N)**2*$dd[$i][$j]/$rr[$i][$j]
			       	- 2*($N_r/$N)*$dr[$i][$j]/$rr[$i][$j];
				$error = ($N/$N_r)**2 * (1 + $ksi) / 
					sqrt($rr[$i][$j]);
			}
			else {
				carp "RR($bins->[$i][0]) or ", 
				"DD($bins->[$i][0]) was 0";
				$ksi = $error = 0;
				
			}
			push @result, [@{$bins->[$i]}, @{$trbins->[$j]}, $ksi, 
					$error];
		}
	}

	return @result;
}
1;

__END__


=head1 NAME

AMIGA::Correlation - Implementation of Landy-Szalay two-point correlation 
estimator.

=head1 SYNOPSIS

	use AMIGA::Correlation qw(calculate_correlation);

=head1 FUNCTIONS 

=head2 calculate_correlation ( ARGS )

Arguments: C<ARGS> -- A hash table with the following possible keys:

=over 4

=item low_bound

Low bound of the spatial distribution of the data. Should be 
an Tools::Vector object.

=item high_bound

High bound of the spatial distribution of the data. Should be 
an Tools::Vector object.

=item data

Data to calculate two-point correlation for. Should be a reference to an
array of Tools::Vector objects.

=item bins 

An array reference of scale bins. Each bin in itself should be an array
reference of two scalars, bin low and high bound. If L<trbins> is given,
then these bins are taken as the parallel component bins.

=item trbins

An array reference of transverse scale bins. If this key is given, then
the correlation is calculated separately for parallel and transverse components.
Parallel scales are given in L<bins>.

=item field

I<Optional.> A precalculated comparison field with the same low and high
bounds as the data, but randomly distributed. Should be a reference to an
array of Tools::Vector objects.

If this key is not specified, then Perl C<rand()> will be used to calculate
a comparison field. This might not be as reliable as using some proven
algorithm, such as the Mersenne Twister.

=item external_prg

I<Optional.> Path and filename of an external program to do calculationally
intensive operations, i.e. pair length calculation and binning, and calculating
the estimator value. Using the provided B<correlator> program is strongly
recommended, as calculating the correlation with pure Perl is rather slow.

The program should conform to the following calling convention:

	external_prg OUT DATA CMP

These command line arguments should be of the following form:

=over 4

=item *

C<OUT> -- Name of the in which to output the data. Output should have
one row per scale bin and each row should have the following values in
order: bin low bound, bin high bound, value of correlation estimator, 
Poisson error.

=item *

C<DATA> -- Name of the data input file. It should have the following rows
in order:

=over 4

=item *

number of scale bins

=item *

scale bin low bound and high bound separated by whitespace for each 
scale bin

=item *

number of data points

=item *

data low bound (3 components separated with whitespace)

=item *

data high bound (3 components separated with whitespace)

=item *

data points (3 components separated with whitespace)

=back

Example of a data file:

	6
	0.01 0.03162
	0.03162 0.1
	0.1 0.3162
	0.3162 1.0
	1.0 3.162
	3.162 10.0
	200
	0 0 0
	1 1 1
	0.290558803201751 0.970668561381416 0.485018308778432
	0.83400878482184 0.216630258260086 0.542955000637434
	.
	.
	etc.

=item * 

C<CMP> -- Name of the comparison field file. It should have the following rows
in order:

=over 4

=item *

number of data points

=item *

data low bound (3 components separated with whitespace)

=item *

data high bound (3 components separated with whitespace)

=item *

data points (3 components separated with whitespace)

=back

Example of a comparison field file:

	1000
	0 0 0
	1 1 1
	0.290558803201751 0.970668561381416 0.485018308778432
	0.83400878482184 0.216630258260086 0.542955000637434
	.
	.

=back

=back

If the correlation is done for both components, then the C<DATA> file should
have two sequential sets of bin rows, first for the parallel, then for the
transverse component.

Return value: A 2D array containing one row for each scale bin, or if in
component mode, C<P*T> rows, where C<P> is the number or parallel scale bins and
C<T> the number of transverse scale bins.

=head1 AUTHOR

Written in 2006 for the benefit of Tuorla Observatory Cosmological Simulations
research group by Pauli Pihajoki.

=head1 COPYRIGHT

Copyright (c) 2006 Pauli Pihajoki. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms
as Perl itself.

=cut

