#!/usr/bin/perl
use strict;
use warnings;

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
    if(/^\s*geom_\S+\s*$/)
    {
	#print $_;
	#exit();
	my $line = <IN>;
	chomp($line);
	$line =~s/^\s+//;
	my $num = 1.0*$line;
	#print "$num\n";
	#exit();
	for(my $i=0;$i<$num;$i++)
	{
	    $line = <IN>;
	    chomp($line);
	    #print "$line\n";
	    #exit();
	    if($line=~m/^\s*(\S+)\s+(\S+)\s+(\S+)/)
	    {
		print "$1 $2 $3\n";
		#exit();
	    }
	}
    }
}
close IN;
