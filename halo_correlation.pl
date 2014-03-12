#!/usr/bin/perl


use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use POSIX qw(log10);

# Do some rather ugly stuff so that we don't have to place the libraries in some
# @INC directory and can just keep all the script files in one place
my (undef, $modpath, undef) = File::Spec->splitpath($0);
push @INC, File::Spec->canonpath($modpath);
require AMIGA::HaloAnalyzer; import AMIGA::HaloAnalyzer 
qw(parse_datafile get_mass_interval get_parameter_interval);

my $outfile = '';
my $box = 10;
my $mass_lb = 0;
my $mass_hb = 0;
my $scale_lb = 1.0;
my $scale_hb = 10.0;
my $bins = 10;
my $log = 1;
my $components = 0;
my $ages = 0;
my $fieldfile = '';

# Issue help message if no arguments.
&usage unless (@ARGV);

GetOptions(
	'out=s' => \$outfile,
	'minmass=f' => \$mass_lb,
	'maxmass=f' => \$mass_hb,
	'box=f' => \$box,
	'field=s' => \$fieldfile,
	'minscale=f' => \$scale_lb,
	'maxscale=f' => \$scale_hb,
	'bins=i' => \$bins,
	'log!' => \$log,
	'components!' => \$components,
	'ages!' => \$ages,
	'help' => \&usage
);

sub usage {
	print "$0 OPTIONS\n";
	print <<HALT;
Valid options:
--out           Output file name
--field         File name of the random sample file
--minmass       Mass low bound 
--maxmass       Mass high bound (default: include halos of any mass)
--box           Box size (in MPc)
--minscale      Scale low bound
--maxscale      Scale high bound
--bins          Number of scale bins to use
--log/--nolog   Whether to have logarithmic scale bins (default: on)
--ages/--noage  Whether to calculate correlation for youngest and oldest 20 % 
                (default: off)
--components    Whether to calculate parallel and transverse correlations 
                separately
HALT
	exit 1;
}

# Check what we got.
die "Invalid output file $outfile" unless ($outfile);
die "Invalid data box $box" if ($box <= 0);
die "Invalid mass limits" if ($mass_lb > $mass_hb);

my @binfield = ();

# If a precalculated binomial field was specified, try to read it in.
if ($fieldfile) {
	# Parse it in, and discard the first 3 rows that contain the
	# preamble.
	print "Reading a binomial random field from $fieldfile...\n";
	@binfield = &parse_datafile($fieldfile);
	@binfield = @binfield[3..$#binfield];
	unless (@binfield) {
		die "Reading data from $fieldfile failed";
	}
}

# Initialize library
my $analyzer = AMIGA::HaloAnalyzer->new(
	verbose => 1, 
	data_path => File::Spec->curdir(), 
	binary_path => $modpath,
	box_size => $box
);

# Construct scale bins according to given commandline options.
my @scalebins;
my ($low, $high);
my $inc;
if ($log) {
	$inc = (log10($scale_hb) - log10($scale_lb))/$bins;
}
else {
	$inc = ($scale_hb - $scale_lb)/$bins;
}
for (my $i=0; $i < $bins; $i++) {
	if ($log) {
		$low = log10($scale_lb) + $i*$inc;
		$high = log10($scale_lb) + ($i+1)*$inc;
		push @scalebins, [10**$low, 10**$high];
	}
	else {
		$low = $scale_lb + $i*$inc;
		$high = $scale_lb + ($i+1)*$inc;
		push @scalebins, [$low, $high];
	}
}

# Get the halo data.
my @halodata = $analyzer->get_all_halodata($analyzer->{z_values}[0]);

# Find the high and low mass indexes, if a mass interval was specified.
# Otherwise get all halos.
my ($hi, $li);
if ($mass_lb == 0 && $mass_hb == 0) {
	$hi = 0;
	$li = $#halodata;
}
else {
	(undef, $hi, $li) = &get_mass_interval(\@halodata, $mass_lb, $mass_hb);
}

# If we want to separate the halos by the formation age, then get halohistory
# for those halos
if ($ages) {
	print "Fetching halohistory for indexes $hi - $li\n";
	my %halohistory = $analyzer->halohistory($hi..$li);

	# Go through the halos file line by line, and add the formation time
	# column
	for (my $i=$hi; $i < $li+1; $i++) {
		# Get formation time for this halo 
		print "Finding formation time for halo $i... ";
		my $age = ($analyzer->find_age($halohistory{$i}))[0];
		print "z = $age\n";
		push @{$halodata[$i]}, $age;
	}
}

# Then crop away the rest of the halodata
@halodata = @halodata[$hi..$li];

my @young;
my @old;
if ($ages) {
	# Then, if necessary, sort ascending by formation redshift and find the
	# youngest and oldest 20 %
	my $num_halos = @halodata;
	my $zi = scalar(@{$halodata[0]})-1;
	@halodata = sort {$a->[$zi] <=> $b->[$zi]} @halodata;
	@young = @halodata[0..($num_halos/5)];
	@old = @halodata[($#halodata-$num_halos/5)..$#halodata];
	printf "Youngest 20%%: %d oldest 20%%: %d\n", scalar(@young), 
		scalar(@old);
}

# Get the correlations for both sets
my @yc;
my @oc;
my %args;
$args{bins} = \@scalebins;
$args{field} = \@binfield if (@binfield);
$args{trbins} = \@scalebins if ($components);
if ($ages) {
	@yc = $analyzer->correlation( halodata => \@young, %args);
	@oc = $analyzer->correlation( halodata => \@old, %args);
}
else {
	@yc = $analyzer->correlation( halodata => \@halodata, %args);
}

# Print it out to the file
open (OUT, ">$outfile") or die "Couldn't open file for reading: $outfile";

my $print_result = sub {
	my $fh = shift;
	my $row = shift;

	my $str = "%-16g " x scalar(@$row);
	$str .= "\n";
	$str = sprintf($str, @$row);
	print $fh $str;
};

if ($ages) {
	print OUT "# Mass interval: [$mass_lb, $mass_hb]\n";
	print OUT "# Youngest 20%:\n";
	print OUT "# r_low r_high ", 
	"ksi([r_low, r_high]) error\n";
	foreach my $row (@yc) {
		&$print_result(*OUT, $row);
	}
	print OUT "\n\n";

	print OUT "# Oldest 20%:\n";
	print OUT "# r_low r_high ", 
	"ksi([r_low, r_high]) error\n";
	foreach my $row (@oc) {
		&$print_result(*OUT, $row);
	}
	print OUT "\n\n";
}
else {
	print OUT "# Mass interval: [$mass_lb, $mass_hb]\n";
	print OUT "# r_low r_high ", 
	"ksi([r_low, r_high]) error\n";
	foreach my $row (@yc) {
		&$print_result(*OUT, $row);
	}
	print OUT "\n\n";
}

close OUT;

