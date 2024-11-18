#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <profile file name with t dpp tuples> <hydrostatic gradient> <top: x y z> <bottom: x y z> opt:<perturbation interval top> <perturbation interval bottom> <scale dpp> <offset t>" if($narg < 8 || $narg==9 || $narg > 13);

my $filename = shift;
die "cannot open $filename" if(!open(IN,"<$filename"));

my $ppgrad = shift;
my @x0 = (0,0,0);
my @x1 = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
    $x1[$i] = shift;
}
for(my $i = 0; $i < 3; $i++)
{
    $x0[$i] = shift;
}

my $zpp1 = $x1[2];
my $zpp0 = $x0[2];

my $nz = 4;
if($narg > 8)
{
    $zpp1 = shift;
    $zpp0 = shift;
    $nz += 2;
}
my $dppfct = 1.0;
my $toffset = 0.0;
if($narg > 10)
{
    $dppfct = shift;
    if($narg > 11)
    {
	$toffset = shift;
    }
}
print "USING: x0 ($x0[0] $x0[1] $x0[2]) 
              x1 ($x1[0] $x1[1] $x1[2]) 
              zpp ($zpp0 $zpp1) 
              dpp ($dppfct) toffset ($toffset)
";

##################################
#PRINT THE COORDINATES IN X AND Y
##################################
{
    open(XOUT,">x");
    open(YOUT,">y");
    my $dx = ($x1[0] - $x0[0]) / 3;
    my $dy = ($x1[1] - $x0[1]) / 3;
    for(my $i = 0; $i < 4; $i++)
    {
	my $xx = $x0[0] + $dx * $i;
	my $yy = $x0[1] + $dy * $i;
	print XOUT "$xx\n";
	print YOUT "$yy\n";
    }
    close XOUT;
    close YOUT;
}


##################################
#PRINT THE COORDINATES IN Z
##################################
my @zs = ();
{
    if($nz == 4)
    {
	open(ZOUT,">z");
	my $dz = ($x1[2] - $x0[2]) / 3;
	for(my $i = 0; $i < 4; $i++)
	{
	    my $zz = $x0[2] + $dz * $i;
	    push(@zs, $zz);
	    print ZOUT "$zz\n";
	}
	close ZOUT;
    }
    else
    {
	my $dz = $zpp1 - $zpp0;
	my $dzsmall = 0.01 * $dz; 
	push(@zs, $x0[2]);
	my $zpp0pre = $zpp0 - $dzsmall;
	push(@zs, $zpp0pre);
	push(@zs, $zpp0);
	push(@zs, $zpp1);
	my $zpp1post = $zpp1 + $dzsmall;
	push(@zs, $zpp1post);
	push(@zs, $x1[2]);
	open(ZOUT,">z");
	for(my $i = 0; $i <= $#zs; $i++)
	{
	    print ZOUT "$zs[$i]\n";
	}
	close ZOUT;
    }
}

##################################
#PRINT THE PORE PRESSURES INTO VOXEL
##################################
open(T,">t");
my @dpps = ();
while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line=~s/^\s+//;
    my @arr  = split(/\s+/,$line);
    my $tt = $arr[0];
    if($tt > $toffset)
    {
	$tt -= $toffset;
    }
    print T "$tt\n";
    my $dpp = $dppfct * $arr[1];
    push(@dpps, $dpp);
}
close IN;
close T;

my $nt = $#dpps + 1;
open(V, ">voxel");
print V "4 4 $nz $nt\n";
for(my $ipp = 0; $ipp <= $#dpps; $ipp++)
{
    for(my $k = 0; $k <= $#zs; $k++)
    {
	my $zz = $zs[$k];
	my $pp = $ppgrad * $zz;
	if($zz >= $zpp0 && $zz <= $zpp1)
	{
	    $pp += $dpps[$ipp];
	    print "   PP: $pp ($zz x $ppgrad + $dpps[$ipp])\n";
	}
	else
	{
	    print "   PP: $pp ($zz x $ppgrad)\n";
	}
	for(my $i = 0; $i < 16; $i++)
	{
	    print V "$pp ";
	}
	print V "\n";
    }
}
close V;
