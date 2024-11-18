#!/bin/perl
use strict;
my $narg = $#ARGV + 1;
die "usage: $0 <filename>\n" if($narg != 1);
my $filename = shift;
die "cannot open $filename\n" if(!open(IN,"<$filename"));
while(<IN>)
{
    if(/^source/)
    {
	my $t = -1;
	if(/t0\=(\S+)/)
	{
	    $t = $1;
	}
	else
	{
	    die "could not find time in entry: $_";
	}
	my $moment = -1;
	if(/m0\=(\S+)/)
	{
	    $moment = (log($1)/log(10) - 9.1)/1.5;
	}
	else
	{
	    die "could not find time in entry: $_";
	}
	print "$t $moment\n";
    }
}
