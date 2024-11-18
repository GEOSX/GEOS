#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <filename>\n" if($narg!=1);

#dimension
#n
#point coordinates
my @pts = ();
my $ff = shift;
open(IN,"<$ff");
while(<IN>)
{
    my $line = $_;
    $line =~s/^\s+//;
    chomp($line);
    my @arr = split(/\s+/,$line);
    push(@pts,\@arr);
}
close IN;

#print the formatted data
my $nn = $#pts + 1;
print "3\n$nn\n";
for(my $i = 0; $i <= $#pts; $i++)
{
    for(my $j = 0; $j < 3; $j++)
    {
	print " $pts[$i][$j]";
    }
    print "\n";
}
