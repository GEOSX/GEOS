#!/usr/bin/perl
use strict;
my $narg = $#ARGV + 1;
if($narg != 4)
{
    die "usage: $0 <filename> <dx> <dy> <dz>\n";
}

my $filename = shift;
my @dx = (0,0,0);
$dx[0] = shift;
$dx[1] = shift;
$dx[2] = shift;

open(IN,"<$filename");
while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line =~s/^\s+//;
    my @arr = split(/\s+/,$line);
    my @x0 = ($arr[0],$arr[1],$arr[2]);
    my @x1 = ($arr[3],$arr[4],$arr[5]);
    my @var = ();
    for(my $i = 0; $i < 3; $i++)
    {
	$x0[$i] += $dx[$i];
	$x1[$i] += $dx[$i];
	print "$x0[$i] ";
    }
    for(my $i = 0; $i < 3; $i++)
    {
	print "$x1[$i] ";
    }

    for(my $i = 6; $i <= $#arr; $i++)
    {
	print "$arr[$i] ";
    }
    print "\n";
}
close IN;
