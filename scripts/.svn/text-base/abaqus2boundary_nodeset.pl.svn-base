#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <abaqus file> <new nodeset name> <x> <y> <z> <nx> <ny> <nz>\nThis script extracts a nodeset from a set of nodes based on a half-space\n" if($narg != 8);

################
# GET ARGUMENTS
################

my $file = shift;
my $newName = shift;
my @x = (0,0,0);
my @nx = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
    $x[$i] = shift;
}
for(my $i = 0; $i < 3; $i++)
{
    $nx[$i] = shift;
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
	    my $dot = 0;
	    for(my $i = 0; $i < 3; $i++)
	    {
		$dot += ($arr[$i+1] - $x[$i]) * $nx[$i];
	    }
	    if($dot >= 0)
	    {
		#print "$arr[1] $arr[2] $arr[3]\n"
		print "    $arr[0],\n";
	    }
	}
    }
}
close IN;
