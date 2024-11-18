#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <fracflow mesh filename>\n" if($narg != 1);

my $fname = shift;

my $ncurr = 1;

print "*HEADING
*NODE,NSET=ALLNODES
";
die "cannot open $fname\n" if ( !open( IN, "<$fname" ) );
while(<IN>)
{
    if(/^\s*nod/)
    {
	my $line = $_;
	$line =~s/^\s+//;
	chomp($line);
	my @arr = split(/\s+/,$line);
	print "$ncurr,$arr[1],$arr[2],0.0\n";
	$ncurr++;
    }
}
close IN;

print "*ELEMENT,TYPE=STRI,ELSET=EB1
";
die "cannot open $fname\n" if ( !open( IN, "<$fname" ) );
my $elcurr = 1;
while(<IN>)
{
    if(/^\s*tri/)
    {
	my $line = $_;
	$line =~s/^\s+//;
	chomp($line);
	my @arr = split(/\s+/,$line);
	print "$elcurr,$arr[1],$arr[2],$arr[3]\n";
	$elcurr++;
    }
}
close IN;
