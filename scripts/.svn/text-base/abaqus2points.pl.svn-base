#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "This script extracts carriage-return delimited tuples consisting of space-delimited nodal positions (3-vectors)
usage: $0 <ABAQUS input filename>
" if($narg != 1);

my $file = shift;
die "$0 : cannot open $file!\n" if(!open(IN,"<$file"));
my $node_on = 0;
while(<IN>)
{
    if(/^\*(\S+)/)
    {
	if(/NODE/)
	{
	    $node_on = 1;
	}
	else
	{
	    $node_on=0;
	}
    }
    else
    {
	if($node_on)
	{
	    my $line = $_;
	    chomp($line);
	    my @arr = split(/\,\s*/,$line);
	    print "$arr[1] $arr[2] $arr[3]\n";
	}
    }
}
close IN;
