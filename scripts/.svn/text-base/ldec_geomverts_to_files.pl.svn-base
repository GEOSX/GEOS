#!/usr/bin/perl
use strict;
use warnings;

my $narg = $#ARGV + 1;
die "usage: $0 <filename> <file_prefix>\n" if($narg != 2);

my $filename = shift;
my $prefix = shift;

open(IN,"<$filename");
my $n = 0;
my $on = 0;
while(<IN>)
{
    if(/^\s*begin_block/)
    {
#geom_\S+\s*$/)
	if($on == 1)
	{
	    close OUT;
	}
	my $str = sprintf(">$prefix%06d",$n);
	$n++;
	open(OUT,"$str");
	$on = 1;
    }
    if($on == 1)
    {
	print OUT $_;
    }
}
if($on == 1)
{
    close OUT;
}
close IN;
