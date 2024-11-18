#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <point pair file>\nfile should be in the following format: x0 y0 z0 x1 y1 z1 scalars ... \n" if($narg != 1);

my $filename = shift;

print "# vtk DataFile Version 1.0\ndata\nASCII\nDATASET POLYDATA\nPOINTS ";
my $num =0;
open(IN,"<$filename");
my $nvars = 0;
while(<IN>)
{
    $num++;
    my $line = $_;
    chomp($line);
    my @arr = split(/\s+/,$line);
    if($nvars <= ($#arr-6))
    {
	$nvars = $#arr - 5;
    }
}
close IN;
$num *= 2;
print "$num double\n";

open(IN,"<$filename");
my @vars = ();
while(<IN>)
{
    my $line = $_;
    chomp($line);
    my @arr = split(/\s+/,$line);
    if($#arr<=4)
    {
	die "insufficient entries!";
    }
    print "$arr[0] $arr[1] $arr[2]\n$arr[3] $arr[4] $arr[5]\n";
    my @extra = ();
    for(my $i = 6; $i <= $#arr; $i++)
    {
	push(@extra, $arr[$i]);
    }
    for(my $i = ($#arr + 1); $i < $nvars; $i++)
    {
	push(@extra,0);
    }
    push(@vars,[@extra]);
}
close IN;

my $numLines = $num/2;
my $size = $numLines * 3;
print "LINES $numLines $size\n";
for(my $i = 0; $i < $numLines; $i++)
{
    my $i0 = 2*$i;
    my $i1 = $i0 + 1;
    print "2 $i0 $i1\n";
}

if($nvars > 0)
{
    print "POINT_DATA $num\n";
    for(my $ii = 0; $ii < $nvars; $ii++)
    {
	my $name = "scalar$ii";
	print "SCALARS $name double\nLOOKUP_TABLE default\n";
	for(my $i = 0; $i < $numLines; $i++)
	{
	    my $val = $vars[$i][$ii];
	    print "$val\n$val\n";
	}
    }
}
if($nvars == 0)
{
    print "POINT_DATA $num\nSCALARS default double\nLOOKUP_TABLE default\n";
    for(my $i = 0; $i < $numLines; $i++)
    {
	print "0\n0\n";
    }
}

exit();
print "\n\n-----------------------------\n\n";
for(my $i = 0; $i < $numLines; $i++)
{
    print "LINE $i ";
    for(my $ii = 0; $ii < $nvars; $ii++)
    {
	my $val = $vars[$i][$ii];
	print "$val ";
    }
    print "\n-----------------------------\n";
}
