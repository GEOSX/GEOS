#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <file name> <scale-axis> <scale-r> <opt: x-axis y-axis z-axis> (narg = $narg currently)\n"
  if ( $narg != 3 && $narg != 6);

#--------------------------------
#READ IN ARGUMENTS
#--------------------------------

my $fname = shift;
die "cannot open $fname\n" if ( !open( IN, "<$fname" ) );

my $scale_axis = shift;
my $scale_r = shift;

#die "ascale: $scale_axis rscale: $scale_r\n";

my @axis = (0,0,1);
if($narg > 3)
{
    my $dd = 0;
    for(my $i = 0; $i < 3; $i++)
    {
	$axis[$i] = shift;
	$dd += $axis[$i] * $axis[$i];
    }
    $dd = sqrt($dd);
    if($dd != 0.0 && $dd != 1.0)
    {
	$dd = 1.0/$dd;
    }
    else
    {
	$dd = 1.0;
    }
    for(my $i = 0; $i < 3; $i++)
    {
	$axis[$i] *= $dd;
    }
}

#--------------------------------
#READ ABAQUS
#--------------------------------
{
	my $node_on = 0;
	my $icurr   = 0;
	while (<IN>) {
		if (/^\s*\*\s*(\S+)/) {
			if (/NODE/) {
				$node_on = 1;
			}
			else {
				$node_on = 0;
			}
			print $_;
		}
		else {
			if ($node_on) {
				my $line = $_;
				chomp($line);
				$line =~ s/^\s+//;
				my @arr = split( /\,\s*/, $line );
				die "node_on but this doesn't look like a node: $_" if ( $#arr != 3 );

				my @x = ( $arr[1], $arr[2], $arr[3] );

				my $dd = 0;
				for(my $i = 0; $i < 3; $i++)
				{
				    $dd += $x[$i] * $axis[$i];
				}
				for(my $i = 0; $i < 3; $i++)
				{
				    my $xaxis = $dd * $axis[$i];
				    $x[$i] -= $xaxis;

				    $x[$i] *= $scale_r;
				    $xaxis *= $scale_axis;

				    $x[$i] += $xaxis;
				}
				$icurr++;

				print "$arr[0], $x[0], $x[1], $x[2]\n";
			}
			else {
				print $_;
			}
		}
	}
}
close IN;
