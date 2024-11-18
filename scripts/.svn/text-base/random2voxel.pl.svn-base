#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <top: x y z> <bottom: x y z> <nx ny nz> <range: min max> ($narg arguments instead)\n" if($narg != 11);

my @x0 = (0,0,0);
my @x1 = (0,0,0);
my @dx = (0,0,0);
my @nn = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
    $x1[$i] = shift;
}
for(my $i = 0; $i < 3; $i++)
{
    $x0[$i] = shift;
    die "$0 : $i -> x0 should be less than x1" if($x0[$i] >= $x1[$i]);
}
for(my $i = 0; $i < 3; $i++)
{
	$nn[$i] = shift;
	if($nn[$i] < 4)
	{
		$nn[$i] = 4;
	}
	$dx[$i] = ($x1[$i] - $x0[$i]) / $nn[$i];
	die "$0 : $i -> dx is <= 0 -> $dx[$i]\n" if($dx[$i] <= 0);
}
my $vmin = shift;
my $vmax = shift;

print "       USING: x0 ($x0[0] $x0[1] $x0[2]) 
              x1 ($x1[0] $x1[1] $x1[2]) 
              dx ($dx[0] $dx[1] $dx[2]) 
              n ($nn[0] $nn[1] $nn[2]) 
              min/max: $vmin $vmax
";

##################################
#PRINT THE COORDINATES IN X AND Y
##################################
{
	for(my $ii = 0; $ii < 3; $ii++)
	{
		my $str = ">";
		if($ii == 0)
		{
			$str .= "x";
		}
		elsif($ii == 1)
		{
			$str .= "y";			
		}
		else
		{
			$str .= "z";			
		}
		die "cannot open $str\n" if(!open(OUT,$str));
		for(my $i = 0; $i < $nn[$ii]; $i++)
		{
		    my $xx = $x0[$ii] + $dx[$ii] * $i;
		    print OUT "$xx\n";	    	
		}
		close OUT;
	}
}

##################################
#PRINT THE RANDOM VALUES INTO VOXEL
##################################
die "cannot open voxel\n" if(!open(V, ">voxel"));
print V "$nn[0] $nn[1] $nn[2]\n";
my $dv = $vmax - $vmin;
for(my $iz = 0; $iz < $nn[2]; $iz++)
{
	my $nnn = $nn[0] * $nn[1];
    for(my $k = 0; $k < $nnn; $k++)
    {
    	my $vv = $vmin + rand() * $dv;
    	print V "$vv ";    	
    }
    print V "\n";
}
close V;
