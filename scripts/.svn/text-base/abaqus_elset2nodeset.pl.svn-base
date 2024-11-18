#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <abaqus file> <elset name> <new nodeset name>\nThis script extracts a nodeset from an elset\n" if($narg != 3);

my $file = shift;
my $elsetName = shift;
my $nsetName = shift;

open(IN,"<$file");
my $elset_on = 0;

my %nodes;

while(<IN>)
{
    if(/^\s*\*(\S+)/)
    {
	if(/ELEMENT/)
	{
	    if(/$elsetName/)
	    {
		$elset_on = 1;
	    }
	    else
	    {
		$elset_on = 0;
	    }
	}
	else
	{
	    $elset_on = 0;
	}
    }
    else
    {
	if($elset_on)
	{
	    my $line = $_;
	    $line=~s/^\s+//;
	    chomp($line);
	    my @arr = split(/\,\s*/,$line);
	    for(my $i = 1; $i <= $#arr; $i++)
	    {
		$nodes{$arr[$i]} = 1;
	    }
	}
    }
}
close IN;

my $index = 0;
print "*NSET, NSET=$nsetName
   ";
foreach (keys %nodes)
{
    print " $_";
    if($index > 3)
    {
	print ",\n   ";
	$index = 0;
    }
    else
    {
	print ",";
	$index++;
    }
}
print "\n";
