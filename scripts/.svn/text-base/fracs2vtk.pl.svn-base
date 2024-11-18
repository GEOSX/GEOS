#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
if($narg != 1)
{
    print "usage: $0 <fracgen output>\n";
    exit();
}
my $filename = shift;

print "# vtk DataFile Version 1.0\ndata\nASCII\nDATASET POLYDATA\nPOINTS ";
my $num =0;
open(IN,"<$filename");
while(<IN>)
{
    $num++;
}
$num *= 2;
close IN;
print "$num float\n";

my @apertures = ();
open(IN,"<$filename");
while(<IN>)
{
    if(/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
    {
	my $x = $1;
	my $y = $2;
	my $dx = $3;
	my $dy = $4;
	my $len = $5;
	my $ap = $6;

	my $lhalf = 0.5 * $len;
	my $x0 = $x - $dx * $lhalf; 
	my $x1 = $x + $dx * $lhalf; 
	my $y0 = $y - $dy * $lhalf; 
	my $y1 = $y + $dy * $lhalf; 
	print "$x0 $y0 0\n$x1 $y1 0\n";
	push(@apertures,$ap);
    }
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

print "POINT_DATA $num\nSCALARS apertures float\nLOOKUP_TABLE default\n";
for(my $i = 0; $i < $num; $i++)
{
    my $val = $apertures[$i];
    print "$val\n$val\n";
}
