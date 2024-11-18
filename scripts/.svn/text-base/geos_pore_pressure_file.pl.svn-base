#!/bin/perl
use strict;
my $narg = $#ARGV + 1;
die "usage: $0 <dynamic_viscosity> <permeability> <x_injector> <y_injector> <z_injector> <a2> <b2> <c2> <q-file> <t-file> <x-file> <y-file> <z-file>\n" if($narg != 13);

my $fct = shift;
my $perm = shift;
$fct /= $perm;

#note: setting permeability to 1 and visc to 4*pi/100 = 0.1256637061436 then q will give p (at 1/100 of the radius)

my @x0 = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
    $x0[$i] = shift;
}

my @r2 = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
    $r2[$i] = shift;
}

my @files = ();
for(my $i = 0; $i < 5; $i++)
{
    $files[$i] = shift;
}

my @q = ();
{
    die "cannot open flow file: $files[0]" if(!open(IN,"<$files[0]"));
    {
	my $i = 0;
	while(<IN>)
	{
	    $q[$i] = 1.0 * $_;
	    $i++;
	}
    }
    close IN;
}

my @t = ();
{
    die "cannot open time file: $files[1]" if(!open(IN,"<$files[1]"));
    {
	my $i = 0;
	while(<IN>)
	{
	    $t[$i] = 1.0 * $_;
	    $i++;
	}
    }
    close IN;
}

die "time file must have the same number of entries as flow rate file" if($#t != $#q);

my @x = ();
{
    die "cannot open x-file: $files[2]" if(!open(IN,"<$files[2]"));
    {
	my $i = 0;
	while(<IN>)
	{
	    $x[$i] = 1.0 * $_;
	    $i++;
	}
    }
    close IN;
}
my @y = ();
{
    die "cannot open y-file: $files[3]" if(!open(IN,"<$files[3]"));
    {
	my $i = 0;
	while(<IN>)
	{
	    $y[$i] = 1.0 * $_;
	    $i++;
	}
    }
    close IN;
}
my @z = ();
{
    die "cannot open z-file: $files[4]" if(!open(IN,"<$files[4]"));
    {
	my $i = 0;
	while(<IN>)
	{
	    $z[$i] = 1.0 * $_;
	    $i++;
	}
    }
    close IN;
}

my $nt = $#t + 1;
my $nx = $#x + 1;
my $ny = $#y + 1;
my $nz = $#z + 1;

print "$nx $ny $nz $nt\n";
for(my $tt = 0; $tt < $nt; $tt++)
{
    my $qc = $q[$tt];
    for(my $kk = 0; $kk < $nz; $kk++)
    {
	my $z2 = $z[$kk]-$x0[2];
	$z2 *= $z2 / $r2[2];
	for(my $jj = 0; $jj < $ny; $jj++)
	{
	    my $y2 = $y[$jj] - $x0[1];
	    $y2 *= $y2 / $r2[1];
	    for(my $ii = 0; $ii < $nx; $ii++)
	    {
		my $x2 = $x[$ii] - $x0[0];
		$x2 *= $x2 / $r2[0];
		my $rc = sqrt($x2 + $y2 + $z2);
		if($rc < 0.01)
		{
		    $rc = 0.01;
		}
		my $ir = 1.0 / (4.0 * 3.14159265359 * $rc);
		
		#get pore pressure elevation from injector
		my $p = $qc * $ir * $fct;
		
		#add hydrostatic
		if($z[$kk] < 0)
		{
		    #$p -= 0.001 * 9.81 * $z[$kk];
		}

		#convert to MPa
		$p *= 1e6;
		print "$p\n";
	    }
	}
    }
}

