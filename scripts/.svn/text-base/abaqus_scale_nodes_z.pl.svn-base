#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <zmin> <zmax> <file with aperture for each point> <mated flag> <file to read>\n" if($narg != 5);

my $zmin = shift;
my $zmax = shift;
my $zmid = 0.5 * ($zmax + $zmin);

my @apertures = ();
{
    my $aps = shift;
    die "cannot open $aps!\n" if(!open(IN, "<$aps"));
    while(<IN>)
    {
	my $line = $_;
	chomp($line);
	my @arr = split(/\s+/,$line);
	push(@apertures, $arr[3]);
    }
    close IN;
}

my $nmid = ($#apertures + 1) / 2;
my $dz = $zmid - $zmin;
my $dir = 0;
{
    my $mated = shift;
    if($mated > 0)
    {
	$dir = 1;
    }
    else
    {
	$dir = -1;
    }
}

{
    my $abaqus = shift;
    die "cannot open $abaqus!\n" if(!open(IN,"<$abaqus"));
    my $node_on = 0;
    my $icurr = 0;
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
	    print $_;
	}
	else
	{
	    if($node_on)
	    {
		my $line = $_;
		chomp($line);
		my @arr = split(/\,\s+/,$line);
		my $z = $arr[3];
		my $strain = $apertures[$icurr]/$dz;
		if($strain < 0)
		{
		    $strain  = 0;
		}
	       
		if($z > $zmid)
		{
		    #upper mesh
		    $strain = 1 - $strain;
		    $z = $strain * ($z - $zmax) + $zmax;
		}
		elsif($z==$zmid)
		{
		    if($icurr >= $nmid)
		    {
			#lower mesh
			$strain = 1 + $dir*$strain;
			$z = $strain * ($z - $zmin) + $zmin;
		    }
		    else
		    {
			#upper mesh
			$strain = 1 - $strain;
			$z = $strain * ($z - $zmax) + $zmax;
		    }
		}
		else
		{
		    #lower mesh
		    $strain = 1 + $dir*$strain;
		    $z = $strain * ($z - $zmin) + $zmin;
		}
		print "$arr[0], $arr[1], $arr[2], $z\n";
		$icurr++;
	    }
	    else
	    {
		print $_;
	    }
	}
    }
    close IN;
}
