#!/bin/perl
use strict;
use Math::Trig;

my $narg = $#ARGV + 1;
die "usage: $0 <filename=x,y,z,area,density> <mean_strike range_strike mean_length range_length percentage>+\n" if($narg < 6);

my $filename = shift;
my $narg_left = $narg - 1;
my $ptotal = 0;
my @sets = ();
while($narg_left > 0)
{
	die "cannot parse command line" if($narg_left < 5);
	my @entry = ();
	for(my $i = 0; $i < 5; $i++)
	{
		$entry[$i] = shift;
	}
	$ptotal += $entry[4];
	push(@sets,\@entry);
	$narg_left -= 5;
}

die "cannot open $filename" if(!open(IN,"<$filename"));
my $ifrac = 0;
while(<IN>)
{
	my @arr = &SplitMe($_);
	my @x0 = ($arr[0], $arr[1], $arr[2]);
	my $dl = 0.5*sqrt($arr[3]);
	my $ff = $arr[3]*$arr[4];
	my $nf = int($ff);
	my $x = rand();
	if($x <= ($ff - $nf))
	{
		$nf++;
	}
	for(my $i = 0; $i < $nf; $i++)
	{
		my @xx = ($x0[0],$x0[1],$x0[2]);
		$xx[0] += $dl*(rand()-0.5);
		$xx[1] += $dl*(rand()-0.5);
		$x = rand();
		my $sum = 0;
		for(my $j = 0; $j <= $#sets; $j++)
		{
			$sum += $sets[$j][4] / $ptotal;
			if($x <= $sum)
			{
				my $strike = $sets[$j][0] + 2.0*(rand()-0.5)*$sets[$j][1];
				my $length = $sets[$j][2] + 2.0*(rand()-0.5)*$sets[$j][3];
				&PrintPointPair($strike, $length, \@xx, $ifrac);
				$ifrac++;
				last;
			}
		}
	}
}
close IN;

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

sub PrintPointPair
{
	my $strike = shift;
	my $length = shift;
	my $xptr = shift;
	my $ifrac = shift;
	
	my $fct = atan(1) / 45.0;
	$strike *= $fct;
	my @x = @$xptr;
	my @svec = (sin($strike), cos($strike), 0);
	$length *= 0.5;
	my @x1 = (0,0,0);
	my @x0 = (0,0,0);
	for(my $i = 0; $i < 3; $i++)
	{
		$x0[$i] = $x[$i] - $length*$svec[$i];
		$x1[$i] = $x[$i] + $length*$svec[$i];
	}
	print "$x0[0] $x0[1] $x0[2] $x1[0] $x1[1] $x1[2] $ifrac\n";
}
