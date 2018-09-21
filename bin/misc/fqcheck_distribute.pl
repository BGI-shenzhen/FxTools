#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin";
use PLOT;

my $usage=<<"USAGE";

	Script   : draw distribution of base and quality along reads from .fqcheck file
	Usage    : perl $0 <read1_fqcheck_file> [read2_fqcheck_file] <-o out_file_prefix>
	Exmple   : perl $0 ./1.fqcheck ./2.fqcheck -o ./

USAGE


my $out_prefix;
GetOptions("o:s"=>\$out_prefix);
die $usage unless ($ARGV[0]);
my @fqcheck = @ARGV;
if (! $out_prefix) {
	$out_prefix = $fqcheck[0];
	$out_prefix =~ s/1\.fqcheck//;
	$out_prefix ||= './';
}

print "fqcheck_distribute:\n";
print "input:\t$_\n" foreach (@fqcheck);
print "output:\t${out_prefix}base.png\n";
print "output:\t${out_prefix}qual.png\n";

my ($base, $qual, $last_cycle) = read_fqcheck (@fqcheck);

our $out = $out_prefix.'base';
our %gnuplot_set;
$gnuplot_set{terminal} = 'postscript portrait color size 8, 5;';
$gnuplot_set{grid}     = 'front lc rgb \'gray\';';
$gnuplot_set{out}      = "'$out.ps';";
$gnuplot_set{title}    = "'Base percentage composition along reads';";
$gnuplot_set{xlabel}   = "'Position along reads';";
$gnuplot_set{ylabel}   = "'Percent';";
$gnuplot_set{xrange}   = "[0:$last_cycle->[-1]+1];";
$gnuplot_set{yrange}   = "[-2:52];";
PLOT::plot_base ($base, $last_cycle);

$out = $out_prefix.'qual';
$gnuplot_set{out}      = "'$out.ps';";
$gnuplot_set{title}    = "'Distribution of qualities';";
$gnuplot_set{xlabel}   = "'Position along reads';";
$gnuplot_set{ylabel}   = "'Quality';";
$gnuplot_set{yrange}   = "[-2:42];";
$gnuplot_set{key}      = "off;";
PLOT::plot_qual ($qual);

sub read_fqcheck {
	my @fqcheck = @_;
	my @base;
	my @qual;
	my @last_cycle;
	my $cycle=0;
	for (my $r=0; $r<@fqcheck; $r++) {
		open FQCHECK, "<$fqcheck[$r]" or die "Error: cannot open $fqcheck[$r]";
		while (my $line = <FQCHECK>) {
			next unless ($line =~ /^base/);
			$cycle++;
			my @value = split /\s+/, $line;
			my @tmp1 = ($cycle, (@value[2 ... 6]));
			push @base, \@tmp1;
			for (my $i=7; $i<@value; $i++) {
				my @tmp2 = ($cycle, ($i-7), $value[$i]);
				push @{$qual[$cycle]}, \@tmp2;
			}
		}
		close FQCHECK;
		$last_cycle[$r] = $cycle;
	}
	return \@base, \@qual, \@last_cycle;
}

__END__
