package PLOT;
use strict;
use warnings;
use FindBin qw($Bin);

############################################################

sub read_data_file {
	my ($input_file, @cols) = @_;
	my $data;
	my $N = 0;
	open IN, "<$input_file" or die "Error: cannot open $input_file";
	while (<IN>) {
		chomp;
		$_ or next;
		my @d = (split "\t", $_);
		for (my $i=0; $i<@cols; $i++) {
			$data->[$i][$N] = $d[$cols[$i]];
		}
		$N++;
	}
	close IN;
	return $data, $N;
}


sub plot_cumulative_distribution {
	my ($data, $P) = @_;
	my @sort = sort {$b<=>$a} @{$data};
	my $N = @sort;
	my $tmp = "$::out.tmp";
	open TMP, ">$tmp" or die;
	for (my $i=0; $i<$N; $i++) {
		print TMP  $sort[$i], "\t", ($i+1)/$N*100, "\n";
	}
	close TMP;
	$$P{$_} = $sort[int($_*$N/100)] foreach (keys %{$P});

	my $plot = "reset;\n";
	while (my($k,$v) = each %::gnuplot_set) {
		$plot .= "set $k $v\n" if ($v);
	}
	$::gnuplot_set{xrange} =~ /\[(\d+).*?(\d+)]/;
	my $x = ($2 - $1)*0.75 + $1;
	$plot .= "set label $_ '$_% >=$$P{$_}' at $x, $_;\n" foreach (keys %{$P});
	$plot .= "plot '$tmp' u 1:2 w lp lc rgb 'red' lw 2 notitle, '-' w i lc rgb 'blue' lw 2 notitle;\n";
	$plot .= "$$P{$_} $_\n" foreach (keys %{$P});
	$plot .= "e\n";
	$plot .= "reset;\n";
	system "$Bin/gnuplot <<END\n$plot\nEND\n";
	system "convert $::out.ps $::out.png";
	system "rm $::out.ps $::out.tmp";
	print "\t$$P{$_}" foreach (sort keys %{$P});
}


sub plot_scatter_2in1 {
	my ($x, $y1, $y2, $title1, $title2) = @_;
	if ($title1) { $title1 = "t '$title1' "; }
	else { $title1 = 'notitle'; }
	if ($title2) { $title2 = "t '$title2' "; }
	else { $title2 = 'notitle'; }

	my $tmp = "$::out.tmp";
	open TMP, ">$tmp" or die;
	for (my $i=0; $i<@{$x}; $i++) {
		print TMP  "$$x[$i]\t$$y1[$i]\t$$y2[$i]\n";
	}
	close TMP;

	my $plot = "reset;\n";
	while (my($k,$v) = each %::gnuplot_set) {
		$plot .= "set $k $v\n" if ($v);
	}
	$plot .= "plot '$tmp' u 1:2 w p lc rgb 'red' lw 2 $title1, '' u 1:3 w p lc rgb 'blue' $title2;\n";
	$plot .= "reset;\n";
	system "$Bin/gnuplot <<END\n$plot\nEND\n";
	system "convert $::out.ps $::out.png";
	system "rm $::out.ps $::out.tmp";
}


sub plot_signal {
	my $signal = shift;
	my $tmp = "$::out.tmp";
	open TMP, ">$tmp" or die;
	foreach (@{$signal}) {
		print TMP join "\t", @{$_};
		print TMP "\n";
	}
	close TMP;

	my $plot = "reset;\n";
	while (my($k,$v) = each %::gnuplot_set) {
		$plot .= "set $k $v\n" if ($v);
	}
	$plot .= "plot '$tmp' u 1:2 w lp t 'A',";
	$plot .= " '' u 1:3 w lp t 'C',";
	$plot .= " '' u 1:4 w lp t 'G',";
	$plot .= " '' u 1:5 w lp t 'T'\n";
	$plot .= "reset;\n";
	system "gnuplot <<END\n$plot\nEND\n";
	system "convert $::out.ps $::out.png";
	system "rm $::out.ps $::out.tmp";
}


sub plot_base {
	my ($base, $last_cycle) = @_;
	print @{$last_cycle},"\n";
	my $tmp = "$::out.tmp";
	open TMP, ">$tmp" or die;
	foreach (@{$base}) {
		print TMP join "\t", @{$_};
		print TMP "\n";
	}
	close TMP;

	my $plot = "reset;\n";
	while (my($k,$v) = each %::gnuplot_set) {
		$plot .= "set $k $v\n" if ($v);
	}
	$plot .= "plot '$tmp' u 1:2 w l lw 3 t 'A',";
	$plot .= " '' u 1:3 w l lw 3 t 'C',";
	$plot .= " '' u 1:4 w l lw 3 t 'G',";
	$plot .= " '' u 1:5 w l lw 3 t 'T',";
	$plot .= " '' u 1:6 w l lw 3 t 'N'";
	if ($last_cycle) {
		pop @{$last_cycle};
		if (@{$last_cycle}) {
			$plot .=", '-' w i lc rgb 'blue' lw 3 notitle;\n";
			$plot .= "$_ 50\n" foreach (@{$last_cycle});
			$plot .= "e";
		}
	}
	$plot .= "\n";
	$plot .= "reset;\n";
	#print $plot,"\n";
	#print $plot,"\n";
	system "gnuplot <<END\n$plot\nEND\n";
	system "convert $::out.ps $::out.png";
#	system "rm $::out.ps $::out.tmp";
}


sub plot_qual {
	my $qual = shift;
	my $tmp = "$::out.tmp";
	open TMP, ">$tmp" or die;
	foreach (@{$qual}) {
		foreach (@{$_}) {
			print TMP join "\t", @{$_};
			print TMP "\n";
		}
		print TMP "\n";
	}
	close TMP;
	my $plot = "reset;\n";
	while (my($k,$v) = each %::gnuplot_set) {
		$plot .= "set $k $v\n" if ($v);
	}
	$plot .= "set cbrange [0:100];\n";
#	$plot .= "set palette defined (0 '\#ffffff', 10 '\#8080ff', 50 '\#0000ff', 90 '\#000080', 100 '\#000000');\n";
	$plot .= "set palette defined (0 '\#ffffff', 10 '\#00ff00', 30 '\#ffff00', 50 '\#ff0000', 70 '\#800000', 100 '\#000000');\n";
	$plot .= "unset colorbox;\n";
	$plot .= "plot '$tmp' u 1:2:(\\\$3/10) w image;\n";
	$plot .= "reset;\n";
	print $plot,"\n";
	system "gnuplot <<END\n$plot\nEND\n";
	system "convert $::out.ps $::out.png";
	#system "rm $::out.ps $::out.tmp";
}


sub plot_error {
	my $tmp = $::out;
	my $plot = "reset;\n";
	while (my($k,$v) = each %::gnuplot_set) {
		$plot .= "set $k $v\n" if ($v);
	}
	$plot .= "plot '$tmp' u 1:2 w i;\n";
	$plot .= "reset;\n";
	system "$Bin/gnuplot <<END\n$plot\nEND\n";
	system "convert $::out.ps $::out.png";
	system "rm $::out.ps";
}


sub plot_insert {
	my $tmp = $::out;
	my $plot = "reset;\n";
	while (my($k,$v) = each %::gnuplot_set) {
		$plot .= "set $k $v\n" if ($v);
	}
	$plot .= "plot '$tmp' u 1:2 w l;\n";
	$plot .= "reset;\n";
	system "$Bin/gnuplot <<END\n$plot\nEND\n";
	system "convert $::out.ps $::out.png";
	system "rm $::out.ps";
}

############################################################

1;
__END__
