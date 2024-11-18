#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <ABAQUS input filename>\n" if($narg =! 1);

my $file = shift;
open(IN,"<$file");
my $node_on = 0;
my $element_on = 0;
my @nodes = ();
while(<IN>)
{
    if(/^\*(\S+)/)
    {
	if(/NODE/)
	{
	    $node_on = 1;
	    $element_on = 0;
	}
	elsif(/ELEMENT/)
	{
	    $element_on = 1;
	    $node_on=0;
	}
	else
	{
	    $element_on = 0;
	    $node_on=0;
	}
    }
    else
    {
	if($node_on)
	{
	    my @arr = &SplitMe($_);
	    $nodes[$arr[0]] = [$arr[1],$arr[2],$arr[3]];
	}
	elsif($element_on)
	{
	    my $ret = &Inverted($_, \@nodes);
	    die "return 0\n" if($ret == 0);
	}
    }
}
close IN;
print "Finished successfully!\n";


sub SplitMe
{
    my $line = shift;
    $line =~s/^\s+//;
    chomp($line);
    my @arr = split(/\,\s+/,$line);
    return @arr;
}

sub Inverted
{
    my $line = shift;
    my $ndpt = shift;
    my @nodes = @{$ndpt};

    my @facenodes = &SplitMe($line);

    #------------------------
    #get the middle
    #------------------------
    my @x0 = (0,0,0);

    my $nnds = $#facenodes;

    die "Number of nodes != 8: $line" if($nnds != 8);

    for(my $i = 1; $i <= $nnds; $i++)
    {
	for(my $j = 0; $j < 3; $j++)
	{
	    $x0[$j] += $nodes[$facenodes[$i]][$j];
	}
    }
    for(my $j = 0; $j < 3; $j++)
    {
	$x0[$j] /= $nnds;
    }

    #print "X0: $x0[0] $x0[1] $x0[2]\n";

    #------------------------
    #get the faces
    #------------------------

    my @faces = ();

    push(@faces,[3,0,2]);
    push(@faces,[2,0,1]);

    push(@faces,[0,3,7]);
    push(@faces,[0,7,4]);

    push(@faces,[4,7,6]);
    push(@faces,[4,6,5]);

    push(@faces,[0,4,1]);
    push(@faces,[1,4,5]);

    push(@faces,[1,5,6]);
    push(@faces,[1,6,2]);

    push(@faces,[3,2,7]);
    push(@faces,[7,2,6]);

    my $ret = 1;
    for(my $i = 0; $i <= $#faces; $i++)
    {
	my $ia = $facenodes[$faces[$i][0]+1];
	my $ib = $facenodes[$faces[$i][1]+1];
	my $ic = $facenodes[$faces[$i][2]+1];

	my @nx1 = (0,0,0);
	my @v0 = (0,0,0);
	my @v1 = (0,0,0);
	for(my $j = 0; $j < 3; $j++)
	{
	    my $a = $nodes[$ia][$j];
	    my $b = $nodes[$ib][$j];
	    my $c = $nodes[$ic][$j];

	    $nx1[$j] = ((1.0/3.0)*($a + $b + $c)) - $x0[$j];

	    $v0[$j] = $b - $a;
	    $v1[$j] = $c - $a;
	}
	my @nx2 = ();
	$nx2[0] = $v0[1]*$v1[2] - $v0[2] * $v1[1];
	$nx2[1] = $v0[2]*$v1[0] - $v0[0] * $v1[2];
	$nx2[2] = $v0[0]*$v1[1] - $v0[1] * $v1[0];
	
	my $sum = 0;
	for(my $j = 0; $j < 3; $j++)
	{
	    $sum += $nx1[$j] * $nx2[$j];
	}

	#print "nx1: $nx1[0] $nx1[1] $nx1[2]\n";
	#print "nx2: $nx2[0] $nx2[1] $nx2[2]\n";

	if($sum > 0)
	{
	    $ret = 0;
	    return $ret;
	    #print "$i -> 0\n";
	}
	else
	{
	    #print "$i -> 1\n";
	}
    }

    return $ret;
}
