#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <ABAQUS input filename> <point,radius file> <opt: search distance factor = 0>\n" if($narg != 2 && $narg != 3);

my $fabaqus = shift;
my $fpts = shift;
my $search = 0.0;
if($narg == 3)
{
    my $tmp = shift;
    $search += $tmp;
}

my @pts = ();
die "cannot open $fpts!\n" if(!open(IN,"<$fpts"));
while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line =~s/^\s+//;
    my @arr = split(/\s+/,$line);
    if($#arr == 3)
    {
	$arr[4] = 1.0;
	for(my $i = 0; $i <= $#pts; $i++)
	{
	    my $r2 = $pts[$i][3]*$pts[$i][4];
	    my $sf = &ScaleFactor(\@arr,$pts[$i],$r2);
	    if($sf < $pts[$i][4])
	    {
		$pts[$i][4] = $sf;
		$arr[4] = $sf;
	    }
	}
	push(@pts,\@arr);
    }
}
close IN;


open(OUT,">test.3D");
{
    print OUT "X Y Z sf\n";
    for(my $i = 0; $i <= $#pts; $i++)
    {
	print OUT "$pts[$i][0] $pts[$i][1] $pts[$i][2] $pts[$i][4]\n";
    }
}
close OUT;


die "cannot open $fabaqus!\n" if(!open(IN,"<$fabaqus"));
my $node_on = 0;
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
	    $node_on = 0;
	}
	print $_;
    }
    else
    {
	if($node_on)
	{
	    my $line = $_;
	    chomp($line);
	    my @arr = split(/\,\s*/,$line);
	    my @pt = ($arr[1],$arr[2],$arr[3],0);

	    my $found = 0;
	    for(my $i = 0; $i <= $#pts; $i++)
	    {
		my $sf = &ScaleFactor(\@pt,$pts[$i],$pts[$i][3]);
		if($sf <= (1+$search) && $sf >= (1-$search))
		{
		    print "$arr[0]";
		    for(my $j = 0; $j < 3; $j++)
		    {
			my $x = $pts[$i][4] * ($pt[$j]-$pts[$i][$j]) + $pts[$i][$j];
			print ", $x";
		    }
		    print "\n";
		    $found = 1;
		    last;
		}
	    }
	    if($found == 0)
	    {
		print $_;
		#die "couldn't find it :-(";
	    }
	}
	else
	{
	    print $_;
	}
    }
}
close IN;


sub ScaleFactor
{
    my $sp1 = shift;
    my $sp2 = shift;
    my $r2 = shift;

    my @s1 = @{$sp1};
    my @s2 = @{$sp2};
    
    my $dd = 0;
    for(my $i = 0; $i < 3; $i++)
    {
	my $tmp = $s1[$i] - $s2[$i];
	$tmp *= $tmp;
	$dd += $tmp;
    }
    my $rr = $s1[3] + $r2;
    if($rr == 0)
    {
	die "$s1[0] $s1[1] $s1[2] $s1[3] :: $s2[0] $s2[1] $s2[2] $r2\n";
    }

    $rr *= $rr;
    $dd /= $rr;
    return sqrt($dd);
}
