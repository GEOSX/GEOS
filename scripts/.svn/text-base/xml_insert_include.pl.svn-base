#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <filename>\n" if($narg != 1);
my $filename = shift;
&ReadFile($filename);

sub ReadFile
{
    my $filename = shift;
    die "cannot open $filename\n" if(!open(my $fh,"<$filename"));
    while(<$fh>)
    {
	if(/Include\s+file\s*\=\s*\"(\S+)\"/)
	{
	    $filename = $1;
	    &ReadFile($filename);
	}
	else
	{
	    print $_;
	}
    }
    close IN;
}
