#!/bin/perl
use strict;
my $narg = $#ARGV + 1;
die "usage: <file> <nu : dip> <nv : strike>\n" if($narg != 3);

my $ff = shift;
my $nu = shift;
my $nv = shift;

open(IN,"<$ff");
my @faults = ();
my @fd = ();
my $cdd = 0;

#get segment lengths and totals
while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line =~s/^\s+//;
    my @arr = split(/\s+/,$line);
    push(@faults,\@arr);
    
    #get segment length
    my $dd = 0.;
    for(my $i = 0; $i < 3; $i++)
    {
		my $d = $arr[3+$i] - $arr[$i];
		$d *= $d;
		$dd += $d;
    }
    $dd = sqrt($dd);
    
    #add to cumulative segment length
    $cdd += $dd;
    
    #add to list of cumulative lengths
    push(@fd, $cdd);
}
close IN;

#sub-sample in the strike direction
my $dx = $cdd / $nv;

my $ifault = 0;
my $vlast = 0;
for(my $i = 0; $i < $nv; $i++)
{
	#get the fault in which the current point along v lies
    my $v = $dx * $i;
    while($ifault < $#faults && $v > $fd[$ifault])
    {
		$vlast = $fd[$ifault];
		$ifault++;
		die "out of bounds: $v $fd[$#faults] $ifault\n" if($ifault > $#faults);
    }
    my $vnext = $fd[$ifault];
    my $vv = $vnext - $vlast;
    
    #get the distance along that fault
    $vv = ($v - $vlast) / $vv;

	#interpolate the point
    my $x = $faults[$ifault][0] + ($faults[$ifault][3] - $faults[$ifault][0])*$vv;
    my $y = $faults[$ifault][1] + ($faults[$ifault][4] - $faults[$ifault][1])*$vv;
    my $nx = $faults[$ifault][1] - $faults[$ifault][4];
    my $ny = $faults[$ifault][3] - $faults[$ifault][0];
    my $nn = $nx* $nx + $ny * $ny;
    $nn = 1./sqrt($nn);
    $nx *= $nn;
    $ny *= $nn;
    for(my $j = 0; $j < $nu; $j++)
    {
		my $u = $dx * $j;
		#print x:strike y:strike z:dip nx ny nz
		print "$x $y $u $nx $ny 0 \n";
    }
}
