#!/usr/bin/perl
use strict;
use warnings;

die "perl $0 < insert file> < prefix [0 default] : 1 read sd information > \n" unless (@ARGV >0);


my $out =shift;
my $prefix=0 || shift;


my $total_map_reads = 0;
my $total_single = 0;
my $total_repeat_pair = 0;
my $total_uniq_pair = 0;
my $total_uniq_low_pair = 0;
my $total_uniq_normal_pair = 0;

my %insert;
my $insert;

my $length1=0;
my $length2;

my $median=0;

open (FI,"$out") or die $!;
my $line=<FI>;
$total_map_reads = (split /\s+/,$line)[-1];

#	print "$line\t$total_map_reads\n";

$line=<FI>;
$total_single =(split /\s+/,$line)[-1];
#	print "$total_single\n";
$line=<FI>;

$total_repeat_pair =(split /\s+/,$line)[-1];
$line=<FI>;
$total_uniq_pair = (split /\s+/,$line)[-1];
$line=<FI>;
$total_uniq_low_pair= (split /\s+/,$line)[-1];
$line=<FI>;
$total_uniq_normal_pair = (split /\s+/,$line)[-1];
$line=<FI>;
$median = (split /\s+/,$line)[-1];

$line=<FI>;

my($Lsd,$Rsd) = $line =~ /SD: -(\d+)\/\+(\d+)/;

#	print "$Lsd\t$Rsd\n";


#######swimming in the sky  and flying in the sea #####

close FI;

if($prefix)
{
    my $lane = $out;
    $lane =~ s/\.soap\.insert//;

    print "$lane\t$median\t$Lsd\t$Rsd\n";
}

my $gnuplot = find_gnuplot ();
my $set_plot = set_plot ($gnuplot);

my $plot = "$gnuplot <<END\n";
$plot .= $set_plot;
$plot .= "set output '$out.ps'\n";
$plot .= "plot '$out' u 1:2 w l\n";
$plot .= "END\n";
system "$plot";
system "convert ps:$out.ps $out.pdf";
unlink "$out.ps";

sub set_plot {
    my $gnuplot = shift;
    my $version = `$gnuplot -V`;
    $version or die "Error: cannot execute $gnuplot";
    my $set_plot;
    if ($version =~ /gnuplot 4\.0/) {
        $set_plot  = "set terminal postscript portrait color\n";
        $set_plot .= "set size 0.914, 0.64\n";
        $set_plot .= "set bmargin 10\n";
        $set_plot .= "set grid\n";
        $set_plot .= "set title 'Distribution of insert size'\n";
        $set_plot .= "set label 1 '         Mapped reads: $total_map_reads'\n";
        $set_plot .= "set label 2 '         Single reads: $total_single'\n";
        $set_plot .= "set label 3 '          Repeat pair: $total_repeat_pair'\n";
        $set_plot .= "set label 4 '            Uniq pair: $total_uniq_pair'\n";
        $set_plot .= "set label 5 '   low frequency pair: $total_uniq_low_pair'\n";
        $set_plot .= "set label 6 'normal frequency pair: $total_uniq_normal_pair'\n";
        $set_plot .= "set label 7 '       Peak: $median'\n";
        $set_plot .= "set label 8 '       SD: -$Lsd/+$Rsd'\n";
        $set_plot .= "set label 1 at graph 0.05, -0.20 font \"Mono,12\"\n";
        $set_plot .= "set label 2 at graph 0.05, -0.25 font \"Mono,12\"\n";
        $set_plot .= "set label 3 at graph 0.05, -0.30 font \"Mono,12\"\n";
        $set_plot .= "set label 4 at graph 0.05, -0.35 font \"Mono,12\"\n";
        $set_plot .= "set label 5 at graph 0.05, -0.40 font \"Mono,12\"\n";
        $set_plot .= "set label 6 at graph 0.05, -0.45 font \"Mono,12\"\n";
        $set_plot .= "set label 7 at graph 0.05, -0.50 font \"Mono,12\"\n";
        $set_plot .= "set label 8 at graph 0.05, -0.55 font \"Mono,12\"\n";
        $set_plot .= "set xlabel 'Insert size'\n";
        $set_plot .= "set ylabel '# Pair reads'\n";
#	$set_plot .= "set xrange [0:]\n";
        $set_plot .= "set yrange [0:]\n";
        $set_plot .= "set key off\n";

    }
    elsif ($version =~ /gnuplot 4\.4/ )
    {
        $set_plot  = "set terminal postscript portrait color\n";
#    $set_plot .= "set size 0.8104, 0.55\n";
        $set_plot .= "set size 0.7804, 0.54\n";    
#     $set_plot .= "set size 0.770, 0.53\n";    
        $set_plot .= "set bmargin 10\n";
        $set_plot .= "set grid\n";
        $set_plot .= "set title 'Distribution of insert size'\n";
        $set_plot .= "set label 1 '         Mapped reads: $total_map_reads'\n";
        $set_plot .= "set label 2 '         Single reads: $total_single'\n";
        $set_plot .= "set label 3 '          Repeat pair: $total_repeat_pair'\n";
        $set_plot .= "set label 4 '            Uniq pair: $total_uniq_pair'\n";
        $set_plot .= "set label 5 '   low frequency pair: $total_uniq_low_pair'\n";
        $set_plot .= "set label 6 'normal frequency pair: $total_uniq_normal_pair'\n";
        $set_plot .= "set label 7 '       Peak: $median'\n";
        $set_plot .= "set label 8 '       SD: -$Lsd/+$Rsd'\n";
        $set_plot .= "set label 1 at graph 0.05, -0.22 font \"Mono,11\"\n";
        $set_plot .= "set label 2 at graph 0.05, -0.27 font \"Mono,11\"\n";
        $set_plot .= "set label 3 at graph 0.05, -0.32 font \"Mono,11\"\n";
        $set_plot .= "set label 4 at graph 0.05, -0.37 font \"Mono,11\"\n";
        $set_plot .= "set label 5 at graph 0.05, -0.42 font \"Mono,11\"\n";
        $set_plot .= "set label 6 at graph 0.05, -0.47 font \"Mono,11\"\n";
        $set_plot .= "set label 7 at graph 0.05, -0.52 font \"Mono,11\"\n";
        $set_plot .= "set label 8 at graph 0.05, -0.57 font \"Mono,11\"\n";
        $set_plot .= "set xlabel 'Insert size'\n";
        $set_plot .= "set ylabel '# Pair reads'\n";
#	$set_plot .= "set xrange [0:]\n";
        $set_plot .= "set yrange [0:]\n";
        $set_plot .= "set key off\n";
    
    }
    else {
        $set_plot = "set terminal postscript portrait color size 6.4, 6.4\n";	
        $set_plot .= "set bmargin 10\n";
        $set_plot .= "set grid\n";
        $set_plot .= "set title 'Distribution of insert size'\n";
        $set_plot .= "set label 1 '         Mapped reads: $total_map_reads'\n";
        $set_plot .= "set label 2 '         Single reads: $total_single'\n";
        $set_plot .= "set label 3 '          Repeat pair: $total_repeat_pair'\n";
        $set_plot .= "set label 4 '            Uniq pair: $total_uniq_pair'\n";
        $set_plot .= "set label 5 '   low frequency pair: $total_uniq_low_pair'\n";
        $set_plot .= "set label 6 'normal frequency pair: $total_uniq_normal_pair'\n";
        $set_plot .= "set label 7 '       Peak: $median'\n";
        $set_plot .= "set label 8 '       SD: -$Lsd/+$Rsd'\n";
        $set_plot .= "set label 1 at graph 0.05, -0.20 font \"Mono,12\"\n";
        $set_plot .= "set label 2 at graph 0.05, -0.25 font \"Mono,12\"\n";
        $set_plot .= "set label 3 at graph 0.05, -0.30 font \"Mono,12\"\n";
        $set_plot .= "set label 4 at graph 0.05, -0.35 font \"Mono,12\"\n";
        $set_plot .= "set label 5 at graph 0.05, -0.40 font \"Mono,12\"\n";
        $set_plot .= "set label 6 at graph 0.05, -0.45 font \"Mono,12\"\n";
        $set_plot .= "set label 7 at graph 0.05, -0.50 font \"Mono,12\"\n";
        $set_plot .= "set label 8 at graph 0.05, -0.55 font \"Mono,12\"\n";
        $set_plot .= "set xlabel 'Insert size'\n";
        $set_plot .= "set ylabel '# Pair reads'\n";
#	$set_plot .= "set xrange [0:]\n";
        $set_plot .= "set yrange [0:]\n";
        $set_plot .= "set key off\n";
    }
    return $set_plot;
}

sub find_gnuplot {
#	my $gnuplot = "/usr/bin/gnuplot" ;       
    my  $gnuplot = "/usr/local/bin/gnuplot" ;
    $gnuplot="/opt/blc/genome/biosoft/gnuplot-4.4.0/bin/gnuplot" unless (-e $gnuplot);
    $gnuplot="/share/project005/xuxun/heweiming/bin/gnuplot-4.4.3/bin/gnuplot" unless (-e $gnuplot);
    $gnuplot="/share/project005/xuxun/heweiming/bin/gnuplot-4.0.0/bin/gnuplot" unless (-e $gnuplot);
    $gnuplot = "/usr/bin/gnuplot" unless (-e $gnuplot);

    return $gnuplot;
}
