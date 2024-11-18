#!/usr/bin/perl
my $narg = $#ARGV + 1;
if($narg < 1 || $narg > 2)
{
    print "usage: $0 <filename> <opt: npts>\n";
    exit();
}
my $filename = shift;
my $npts = 2;
if($narg > 1)
{
    $npts = shift;
    if($npts < 2)
    {
	$npts = 2;
    }
}
my $nseg = $npts - 1;

open(IN, "<$filename");
while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line =~s/^\s+//;
    my @arr = split(/\s+/,$line);
    my @x0 = ($arr[0],$arr[1],$arr[2]);
    my @x1 = ($arr[3],$arr[4],$arr[5]);
    my @dx = (0,0,0);
    for(my $i = 0; $i <= $#x1; $i++)
    {
	$dx[$i] = ($x1[$i] - $x0[$i])/$nseg;
    }
    for(my $i = 1; $i <= $nseg; $i++)
    {
	for(my $j = 0; $j <= $#dx; $j++)
	{
	    my $xcurr = $x0[$j] + $dx[$j] * $i;
	    print "$xcurr ";
	}
	if(0)
	{
	for(my $j = 6; $j <= $#arr; $j++)
	{
	    print "$arr[$j] ";
	}
    }
	print "\n";
    }
}
close IN;
