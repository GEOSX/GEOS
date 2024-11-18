#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <abaqus file> <nodeset name> <new nodeset name> <x> <y> <z>\nThis script extracts a nodeset from another nodeset\n" if($narg != 6);

my $file = shift;
my $nsetName = shift;
my $newName = shift;
my @x = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
    $x[$i] = shift;
}

open(IN,"<$file");
my $node_on = 0;
my $nset_on = 0;
my @pts = ();
my @nset = ();

my @xmin = (1e100,1e100,1e100);
my @xmax = (-1e100,-1e100,-1e100);

while(<IN>)
{
    if(/^\s*\*(\S+)/)
    {
	if(/NODE/)
	{
	    $node_on = 1;
	}
	else
	{
	    $node_on = 0;
	}
	if(/NSET/)
	{
	    if(/$nsetName/)
	    {
		$nset_on = 1;
	    }
	    else
	    {
		$nset_on = 0;
	    }
	}
	else
	{
	    $nset_on = 0;
	}
    }
    else
    {
	if($node_on)
	{
	    my $line = $_;
	    $line=~s/^\s+//;
	    chomp($line);
	    my @arr = split(/\,\s*/,$line);
	    $pts[$arr[0]] = [$arr[1],$arr[2],$arr[3]];
	}
	elsif($nset_on)
	{
	    my $line = $_;
	    chomp($line);
	    $line=~s/^\s+//;
	    my @arr = split(/\,\s*/,$line);
	    push(@nset,@arr);
	    for(my $ii = 0; $ii <= $#arr; $ii++)
	    {
		for(my $i = 0; $i < 3; $i++)
		{
		    if($xmin[$i]>$pts[$arr[$ii]][$i])
		    {
			$xmin[$i] = $pts[$arr[$ii]][$i];
		    }
		    if($xmax[$i]<$pts[$arr[$ii]][$i])
		    {
			$xmax[$i] = $pts[$arr[$ii]][$i];
		    }
		}
	    }
	}
    }
}
close IN;

if($nsetName eq "ALLNODES")
{
    @nset = ();
    for(my $ii = 1; $ii <= $#pts; $ii++)
    {
	push(@nset,$ii);
	for(my $i = 0; $i < 3; $i++)
	{
	    if($xmin[$i]>$pts[$ii][$i])
	    {
		$xmin[$i] = $pts[$ii][$i];
	    }
	    if($xmax[$i]<$pts[$ii][$i])
	    {
		$xmax[$i] = $pts[$ii][$i];
	    }
	}
    }
}

print "*NSET, NSET=$newName\n";
#die "$x[0] $x[1] $x[2]\n";
for(my $ii = 0; $ii <= $#nset; $ii++)
{
    my $write_out = 1;
    for(my $i = 0; $i < 3; $i++)
    {
	my $xx = $pts[$nset[$ii]][$i];
	if($x[$i] > 0 && $xx != $xmax[$i]) 
	{
	    #print "xx $i $x[$i] $xx != $xmax[$i]\n";
	    $write_out = 0;
	}
	elsif($x[$i] < 0 && $xx != $xmin[$i])
	{
	    #print "yy $i $x[$i] $xx != $xmin[$i]\n";
	    $write_out = 0;
	}
    }
    if($write_out == 1)
    {
	print "    $nset[$ii],\n";
        #$xmax[2] --> $pts[$nset[$ii]][2],\n";
    }
}
