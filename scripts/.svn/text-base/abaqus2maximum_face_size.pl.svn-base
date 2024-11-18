#!/usr/bin/perl
use strict;
use Math::Trig;

my $narg = $#ARGV + 1;
die "usage: $0 <ABAQUS input filename>\n" if($narg != 1);

my $file = shift;
die "unable to open $file" if(!open(IN,"<$file"));

my $node_on = 0;
my $element_on = 0;
my $endstr = "";
my $end_on = 0;

my @nodes = ();
#my @faces = ([0,3,2,1],[4,5,6,7],[0,4,7,3],[1,2,6,5],[2,3,7,6],[0,1,5,4]);
my $dx_max = 0.0;
while(<IN>)
{
    if($end_on == 1)
    {
	$endstr .= $_;
    }
    elsif(/^\*(\S+)/)
    {
	if(/NODE/)
	{
	    $node_on = 1;
	}
	else
	{
	    $node_on=0;
	}
	if(/PROPERTIES/)
	{
	    $endstr .= $_;
	    $end_on = 1;
	}
	elsif(/ELEMENT/)
	{
	    $element_on = 1;
	}
	else
	{
	    $element_on = 0;
	}
    }
    else
    {
	my $line = $_;
	chomp($line);
	my @arr = split(/\,\s+/,$line);
	
	if($node_on)
	{
	    $nodes[$arr[0]] = [$arr[1],$arr[2],$arr[3]];
	}
	elsif($element_on)
	{
	    my @element = ($arr[1],$arr[2],$arr[3],$arr[4],$arr[5],$arr[6],$arr[7],$arr[8]);
	    
	    for(my $i = 0; $i <= $#element; $i++)
	    {
		for(my $j = ($i+1); $j <= $#element; $j++)
		{
		    my $dx = 0.0;
		    for(my $k = 0; $k < 3; $k++)
		    {
			my $tmp = $nodes[$element[$i]][$k] - $nodes[$element[$j]][$k];
			$tmp *= $tmp;
		    }
		    if($dx > $dx_max)
		    {
			$dx_max = $dx;
		    }
		}
	    }
	}
    }
}
close IN;

$dx_max = sqrt($dx_max);
print $dx_max;
