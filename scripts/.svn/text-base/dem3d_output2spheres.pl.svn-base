#!/usr/bin/perl

my $narg = $#ARGV + 1;
if($narg != 1)
{
    print "usage: $0 <filename>\n";
    exit();
}

my $filename = shift;
open(IN,"<$filename");
while(<IN>)
{
    if(/^SPHERE\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)/)
    {
	print "$1 $2 $3 $4\n";
    }
}
close IN;

