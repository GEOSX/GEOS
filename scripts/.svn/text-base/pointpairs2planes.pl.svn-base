#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
if($narg != 4)
{
    print "usage: $0 <filename> <t0 t1 t2>\n";
    exit();
}
my $filename = shift;
my @t = ();
for(my $i = 0; $i < 3; $i++)
{
    my $v = shift;
    push(@t, $v);
}

open(IN, "<$filename");
while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line =~s/^\s+//;
    my @arr = split(/\s+/,$line);
    my @v = ();

    #get the mid-point
    for(my $i = 0; $i < 3; $i++)
    {
	my $x0 = $arr[$i];
	my $x1 = $arr[$i + 3];
	my $xx = 0.5 * ($x0 + $x1);
	print "$xx ";
	my $vv = $x1 - $x0;
	push(@v, $vv);
    }

    #print "\n#<$v[0] $v[1] $v[2]> x <$t[0] $t[1] $t[2]>\n";
    
    #get v x t = n
    my @n = ();
    {
	my $nsum = 0;
	for(my $i = 0; $i < 3; $i++)
	{
	    my $ii = $i-1;
	    my $jj = $i+1;
	    if($ii < 0)
	    {
		$ii = 2;
	    }
	    elsif($jj > 2)
	    {
		$jj = 0;
	    }
	    my $nn = $v[$ii] * $t[$jj] - $v[$jj] * $t[$ii];
	    push(@n, $nn);
	    $nsum += $nn * $nn;
	    #print "$nn ";
	}
	if($nsum != 0 && $nsum != 1)
	{
	    $nsum = 1.0/sqrt($nsum);
	}
	for(my $i = 0; $i <= $#n; $i++)
	{
	    $n[$i] *= $nsum;
	    print "$n[$i] ";
	}
    }
    print "\n";
}
close IN;
