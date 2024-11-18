#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <abaqus file> <new nodeset name> [<x> <y> <z> <nx> <ny> <nz>]+\nThis script extracts a nodeset from a set of nodes based on a series of half-space clips\n" if($narg < 8);

################
# GET ARGUMENTS
################

my $file = shift;
my $newName = shift;
my @xs = ();
my @nxs = ();
for(my $ii = 0; $ii < ($narg-2); $ii+=6)
{
    my @x = (0,0,0);
    for(my $i = 0; $i < 3; $i++)
    {
	$x[$i] = shift;
    }
    my @nx = (0,0,0);
    for(my $i = 0; $i < 3; $i++)
    {
	$nx[$i] = shift;
    }
    push(@xs, \@x);
    push(@nxs, \@nx);
}


################
# PRINT NODESET
################

print "*NSET, NSET=$newName\n";
die "cannot open $file\n" if(!open(IN,"<$file"));
my $node_on = 0;
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
    }
    else
    {
	if($node_on)
	{
	    my $line = $_;
	    $line=~s/^\s+//;
	    chomp($line);
	    my @arr = split(/\,\s*/,$line);
	    my $ok = 1;
	    for(my $ii = 0; $ii <= $#xs; $ii++)
	    {
		#print "x: $xs[$ii][0] $xs[$ii][1] $xs[$ii][2]\nnx: $nxs[$ii][0] $nxs[$ii][1] $nxs[$ii][2]\n";
		my $dot = 0;
		for(my $i = 0; $i < 3; $i++)
		{
		    $dot += ($arr[$i+1] - $xs[$ii][$i]) * $nxs[$ii][$i];
		}
		if($dot >= 0)
		{
		    $ok = 0;
		}
	    }
	    if($ok == 1)
	    {
		print "    $arr[0],\n";
	    }
	}
    }
}
close IN;
