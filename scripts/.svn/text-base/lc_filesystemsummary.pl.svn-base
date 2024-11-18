#!/usr/bin/perl
use strict;
my $narg = $#ARGV + 1;
die "usage: $0\n" if($narg != 0);

my $on = 0;
die "cannot find /usr/local/docs/FileSystemUse\n" if(!open(IN,"</usr/local/docs/FileSystemUse"));
my $max = 0;
my $lmax = "";
while(<IN>)
{
    if(/^\s*Filesystem/)
    {
	$on = 1;
    }
    elsif($on == 1)
    {
	my $line = $_;
	$line =~s/^\s+//;
	chomp($line);
	my @arr = split(/\s+/,$line);
	if($#arr == 4)
	{
	    print "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\n";
	    if($arr[4]=~m/\d+/)
	    {
		if($arr[4] > $max)
		{
		    $lmax = $arr[0];
		    $max = $arr[4];
		}
	    }
	}
	if(/^SCF/)
	{
	    print "**best_through_put: $lmax\n";
	    exit();
	}
    }
}
close IN;
