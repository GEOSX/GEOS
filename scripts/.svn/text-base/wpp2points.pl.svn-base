#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <WPP source definition filename>\n" if($narg != 1);
my $filename = shift;

die "cannot open $filename" if(!open(IN,"<$filename"));
while(<IN>)
{
    my @arr = split(/\s+/);
    my $i = 2;

    my @key=split(/\=/,$arr[$i]);
    $i++;
    my $t0 = $key[1];

    @key=split(/\=/,$arr[$i]);
    $i++;
    my $x=$key[1];

    @key=split(/\=/,$arr[$i]);
    $i++;
    my $y=$key[1]; 

    @key=split(/\=/,$arr[$i]);
    $i++;
    my $z=$key[1]; 

    @key=split(/\=/,$arr[$i]);
    $i++;
    my $m0=$key[1]; 

    print "$x $y $z $m0 $t0\n";
}
close IN;
