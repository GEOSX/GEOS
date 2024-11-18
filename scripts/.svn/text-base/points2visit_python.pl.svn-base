#!/usr/bin/perl
use strict;
use warnings;

my $narg = $#ARGV + 1;
die "usage: $0 <filename>\n" if($narg != 1);

#get data limits
my $filename = shift;
open(IN,"<$filename");
my $num = 0;
my @minmax = (1e100, 1e100, 1e100, 1e100, -1e100, -1e100, -1e100, -1e100);
my $nfields = 4;
while(<IN>)
{
    my $line = $_;
    $line =~s/^\s+//;
    chomp($line);
    my @arr = split(/\s+/,$line);
    for(my $i = 0; $i < $nfields; $i++)
    {
	if($arr[$i]<$minmax[$i])
	{
	    $minmax[$i] = $arr[$i];
	}
	if($arr[$i]>$minmax[$i+$nfields])
	{
	    $minmax[$i+$nfields] = $arr[$i];
	}
    }
    $num++;
}
close IN;

#get approximate molecule radius
my $r = 1.;
{
    for(my $i = 0; $i < 3; $i++)
    {
	$r *= $minmax[$i+$nfields] - $minmax[$i];
    }
    $r /= $num;
    $r = 0.5*($r**(1./3.));
}

#get approximate visualization paramters
my @ctr = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
    $ctr[$i] = 0.5 * ($minmax[$i+$nfields] + $minmax[$i]);
}
my @nrm = (0,0,1);
my @up = (0,1,0);
my $near = $minmax[2] - 0.1*($minmax[$nfields+2] - $minmax[2]);
my $far = $minmax[$nfields + 2] + 0.1*($minmax[$nfields+2] - $minmax[2]);

#print the simple python script
print "import sys

#***SECTION: RESTORE SESSION
RestoreSession(\"default_molecule.session\", 1) # 1 means in .visit directory
ReplaceDatabase(\"$filename.3D\")

#***SECTION: MOVIE
v0=View3DAttributes()
v0.viewNormal = ($nrm[0],$nrm[1],$nrm[2])
v0.focus = ($ctr[0],$ctr[1],$ctr[2])
v0.viewUp = ($up[0],$up[1],$up[2])
v0.viewAngle = 30
v0.parallelScale = 1
v0.nearPlane = $near
v0.farPlane = $far
v0.imagePan = (0,0)
v0.imageZoom = 1
v0.perspective = 1
v0.eyeAngle = 2
v0.centerOfRotationSet = 1
v0.centerOfRotation = ($ctr[0],$ctr[1],$ctr[2])
SetView3D(v0)

#***SECTION: PLOT ATTRIBUTES
p = MoleculeAttributes()
p.elementColorTable = \"hot\"
p.radiusFixed = $r
SetPlotOptions(p)
#ResetPlotOptions(\"Molecule\")

s = SaveWindowAttributes()
s.width = 300
s.height = 300
s.fileName = \"$filename\"
s.format = s.PNG
DrawPlots()
SetSaveWindowAttributes(s)
SaveWindow()
sys.exit(\"\")
";
