#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <filename> <opt: scalar> <opt: dx>\n" if($narg < 1 || $narg > 3);

my $filename = shift;
my $scalar = 0;
my $dx = 0.0;
my $fdx = 0;
if($narg > 1)
{
    $scalar = shift;
    if($narg > 2)
    {
    	$fdx = 1;
    	$dx = shift;
    }
}

#If we need to, recalculate dx to get an integer fit of the length
if($fdx > 0)
{
	my $dd = 0;
	open(IN,"<$filename");
	my $last_line = "";
	while(<IN>)
	{
	    if(/^\#/)
	    {
	    }
	    else
	    {
			my @arr = &SplitMe($_);
			if($#arr == 1)
			{
				$arr[2] = 0;
			}
			if($last_line ne "")
			{
				my @arr_last = split(/\s+/,$last_line);
				my $cd = &Length(\@arr_last, \@arr);
				$dd += $cd;
			}
			$last_line = "$arr[0] $arr[1] $arr[2]";
	    }
	}
	close IN;
	
	my $nn = int($dd / $dx);
	my $dx2 = $dd / ($nn+1);
	if($dx2 < 0.8 * $dx)
	{
		$dx = $dd / $nn;
	}
	else
	{
		$dx = $dx2;
	}
}

#Now, do the real calculation
{
	open(IN,"<$filename");
	my $last_line = "";
	my $last_line_dx = "";
	while(<IN>)
	{
	    if(/^\#/)
	    {
	    }
	    else
	    {
			my @arr = &SplitMe($_);
			if($#arr == 1)
			{
				$arr[2] = 0;
			}
			my $line = "$arr[0] $arr[1] $arr[2]";
			if($last_line ne "")
			{
				if($fdx)
				{
					$last_line_dx = &SubSampleSegment($last_line, $last_line_dx, $line, $dx, $scalar);				
				}
				else
				{
				    print "$last_line $line $scalar\n";
				}
			}
			$last_line = $line;
	    }
	}
	close IN;
	
	if($last_line_dx ne "")
	{
		my @arr0 = &SplitMe($last_line_dx);
		my @arr1 = &SplitMe($last_line);
		my $dd = &Length(\@arr0, \@arr1);
		if($dd > 0.5 * $dx)
		{
			print "$last_line_dx $last_line $scalar\n";
		}
	}
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

sub SubSampleSegment
{
	my $last_line = shift;
	my $last_line_dx = shift;
	my $line = shift;
	my $dx = shift;
	my $scalar = shift;
	
	my @arr_curr = split(/\s+/, $line);
	my @arr_last = split(/\s+/, $last_line);
	
	#if there is a remainder from last time, get the length remaining
	my $dd_ldx = 0;
	my @arr_dx = ();
	if($last_line_dx ne "")
	{
		@arr_dx = split(/\s+/,$last_line_dx);
		$dd_ldx = &Length(\@arr_last, \@arr_dx);
	}
	
	#get the length of the current segment
	my $dd_c = &Length(\@arr_last, \@arr_curr);
	
	#if the sum of the remainder length and current length is not greater than dx,
	#you will need to wait until there is enough cumulative length to get a point
	#so return the last line
	if(($dd_ldx + $dd_c) < $dx)
	{
		if($last_line_dx ne "")
		{
			return $last_line_dx;
		}
		else
		{
			return $last_line;
		}
	}
	
	#get the terminal point on the current segment
	my $dd_offset = 0;
	if($dd_ldx > 0)
	{
		$dd_offset = $dx - $dd_ldx;
	}
	my @x0 = &InterpolatePoint(\@arr_last, \@arr_curr, $dd_offset/$dd_c);	
	#if you had a remainder from last time, write out the first segment
	if($dd_ldx > 0)
	{
	    print "$arr_dx[0] $arr_dx[1] $arr_dx[2] $x0[0] $x0[1] $x0[2] $scalar\n";		
	}
	
	#get the new length to sub-discretize
	my @x1 = ($arr_curr[0], $arr_curr[1], $arr_curr[2]);
	my $dd_curr = &Length(\@x1, \@x0);
	my $n = int($dd_curr / $dx);
	$last_line_dx = "$x0[0] $x0[1] $x0[2]";
	for(my $ii = 1; $ii <= $n; $ii++)
	{
		my $f0 = ($dx * ($ii - 1)) / $dd_curr;
		my $f1 = ($dx * $ii) / $dd_curr;
		my @xx0 = &InterpolatePoint(\@x0, \@x1, $f0);
		my @xx1 = &InterpolatePoint(\@x0, \@x1, $f1);
		print "$xx0[0] $xx0[1] $xx0[2] $xx1[0] $xx1[1] $xx1[2] $scalar\n";
		$last_line_dx = "$xx1[0] $xx1[1] $xx1[2]";
	}
	return $last_line_dx;
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
