#!/bin/perl
use strict;
my $narg = $#ARGV + 1;
die "usage $0: <file> <origin: x, y, z>\n" if($narg != 4);

my $ff = shift;
my @x0 = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
	$x0[$i] = shift;	
}
#get the largest radial distance from the point
my @pps = ();
die "cannot open $ff\n" if(!open(IN, "<$ff"));
my $line = "";
my @cog = ();
my $ddc = 0;
while(<IN>)
{
	$line = $_;
	my @arr = &SplitMe($line);
	push(@pps,\@arr);
	my @pt0 = ($arr[0], $arr[1], $arr[2]);
	my @pt1 = ($arr[3], $arr[4], $arr[5]);
	my $dd = &Length(\@pt0, \@pt1);
	$ddc += $dd;
	my @pt = ();
	for(my $i = 0; $i < 3; $i++)
	{
		$cog[$i] += 0.5 * ($pt0[$i] + $pt1[$i]) * $dd;
	}
}
close IN;

#get the cog and displacement
$cog[0] /= $ddc;
$cog[1] /= $ddc;
$cog[2] /= $ddc;
my @disp = ();
for(my $j = 0; $j < 3; $j++)
{
	$disp[$j] = $x0[$j] - $cog[$j];
}
#die "$disp[0] $disp[1] $disp[2]\n";

#get the values of the scalar for each point pair and print
for(my $i = 0; $i <= $#pps; $i++)
{
	my $ptr = $pps[$i];
	my @arr = @$ptr;
	for(my $j = 0; $j < 3; $j++)
	{
		$arr[$j] += $disp[$j];
		$arr[$j+3] += $disp[$j];
	}
	print "$arr[0]";
	for(my $j = 1; $j <= $#arr; $j++)
	{
		print " $arr[$j]";
	}
	print "\n";
}

sub SplitMe
{
	my $line = shift;
	chomp($line);
	$line =~s/^\s+//;
	$line =~s/\,/ /;
	$line =~s/\s+/ /;
	my @arr = split(/\s+/,$line);
	return @arr;
}

sub Length
{
	my $ptr0 = shift;
	my $ptr1 = shift;
	my @arr0 = @$ptr0;
	my @arr1 = @$ptr1;
	my $dd = 0;
	for(my $i = 0; $i < 3; $i++)
	{
		my $d = $arr1[$i] - $arr0[$i];
		$d *= $d;
		$dd += $d;
	}
	$dd = sqrt($dd);
	return $dd;
}

sub InterpolatePoint
{
	my $ptr0 = shift;
	my $ptr1 = shift;
	my @arr0 = @$ptr0;
	my @arr1 = @$ptr1;
	my $fct = shift;
	
	my @ret = ();
	for(my $i = 0; $i < 3; $i++)
	{
		$ret[$i] = $arr0[$i] + ($arr1[$i] - $arr0[$i]) * $fct;
	}
	return @ret;
}



