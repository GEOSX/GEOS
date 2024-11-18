#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <filename>\n" if($narg != 1);

my $ff = shift;
open(IN,"<$ff");
my $cdd = 0.;
while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line =~s/^\s+//;
    my @arr = split(/\s+/,$line);
    my $dd = 0.;
    for(my $i = 0; $i < 3; $i++)
    {
	my $d = $arr[$i+3] - $arr[$i];
	$d *= $d;
	$dd += $d;
    }
    $dd = sqrt($dd);
    $cdd += $dd;
    print "$dd $cdd\n";
}
close IN;
