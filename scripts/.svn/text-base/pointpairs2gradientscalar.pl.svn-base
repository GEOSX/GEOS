#!/bin/perl
use strict;
my $narg = $#ARGV + 1;
die "usage $0: <file> <scalar0> <scalar1> <origin: x, y, z>\n" if($narg != 6);

my $ff = shift;
my $v0 = shift;
my $v1 = shift;
my @x0 = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
	$x0[$i] = shift;	
}

#get the largest radial distance from the point
my @pps = ();
my $dd_max = 0;
die "cannot open $ff\n" if(!open(IN, "<$ff"));
my $line = "";
while(<IN>)
{
	$line = $_;
	my @arr = &SplitMe($line);
	push(@pps,\@arr);
	my @pt = ($arr[0], $arr[1], $arr[2]);
	my $dd = &Length(\@pt, \@x0);
	if($dd > $dd_max)
	{
		$dd_max = $dd;
	}
}
close IN;
{
	my @arr = &SplitMe($line);
	my @pt = ($arr[3], $arr[4], $arr[5]);
	my $dd = &Length(\@pt, \@x0);
	if($dd > $dd_max)
	{
		$dd_max = $dd;
	}
}

#get the values of the scalar for each point pair and print
for(my $i = 0; $i <= $#pps; $i++)
{
	my $ptr = $pps[$i];
	my @arr = @$ptr;
	for(my $j = 0; $j <= $#arr; $j++)
	{
		print "$arr[$j] ";
	}
	my @pt = ();
	for(my $j = 0; $j < 3; $j++)
	{
		$pt[$j] = 0.5 * ($arr[$j] + $arr[$j+3]);
	}
	my $dd = &Length(\@pt, \@x0);
	my $vv = $dd / $dd_max;
	$vv *= ($v1 - $v0);
	$vv += $v0;
	print "$vv\n";
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





