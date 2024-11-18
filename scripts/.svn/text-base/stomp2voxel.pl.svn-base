#!/bin/perl
use strict; 

my $narg = $#ARGV + 1;
die "usage: $0 <filename>\n" if($narg != 1);
my $filename = shift;
die "cannot open $filename" if(!open(IN,"<$filename"));

my @n = ();
my @vars = ();

my $vindex = 0;
my $vnameindex = 0;
my $nvals = 0;
my $nn = 0;
my $read = -1;

while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line=~s/^\s+//;
    if(/^I\=(\d+)\,\s+J\=(\d+)\,\s+K\=(\d+)\,\s+/)
    {
	push(@n, $1);
	push(@n, $2);
	push(@n, $3);
	$nn = $n[0] * $n[1] * $n[2];
	$nvals = $nn;
    }
    elsif(/^VARIABLES/)
    {
	my @arr = split(/\"/, $line);
	for(my $i = 1; $i <= $#arr; $i += 2)
	{
	    push(@vars, $arr[$i]);
	}
    }
    elsif(/^ZONE/)
    {
	#do nothing
    }
    else
    {
	my @arr = split(/\s+/,$line);
	my $nnew = $#arr  + 1;
	if($nvals + $nnew > $nn)
	{
	    $nnew += $nvals;
	    for(my $i = 0; $nvals < $nn; $i++, $nvals++)
	    {
		die "$i $nvals < $nn\n";
		#WRITE ME
		&WriteEntry($arr[$i], $read, \@n, $nvals);
	    }
	    my $fname = "";
	    if($vars[$vindex]=~m/Easting/)
	    {
		$fname = "x";
		$read = 0;
	    }
	    elsif($vars[$vindex]=~m/Northing/)
	    {
		$fname = "y";
		$read = 1;
	    }
	    elsif($vars[$vindex]=~m/Elevation/)
	    {
		$fname = "z";
		$read = 2;
	    }
	    else
	    {
		$read = 3;
		$fname = sprintf("variable_%02d", $vnameindex);
		++$vnameindex;
	    }
	    if($vnameindex > 1)
	    {
		print OUT "VARIABLE = \"$vars[$vnameindex+2]\"\n";
	    }
	    ++$vindex;
	    close OUT;
	    die "cannot open $fname for writing" if(!open(OUT, ">$fname"));
	    #if($vindex > 4)
	    #{
	    #die "just wrote file with $nvals lines ... moving on to $fname\n";
	    #}
	    my $nnext = $nnew - $nvals;
	    for(my $i = $0; $nvals < $nnew; $i++, $nvals++)
	    {
		#WRITE ME
		&WriteEntry($arr[$i], $read, \@n, $nvals);
	    }
	    $nvals = $nnext;
	}
	else
	{
	    for(my $i = 0; $i < $nnew; $i++, $nvals++)
	    {
		#WRITE ME
		&WriteEntry($arr[$i], $read, \@n, $nvals);
	    }
	}
    }
}
close IN;
if($vnameindex > 1)
{
    print OUT "VARIABLE = \"$vars[$vnameindex+2]\"\n";
}
close OUT;

sub WriteEntry
{
    my $val = shift;
    my $read = shift;
    my $n_ptr = shift;
    my $nvals = shift;

    my @nn = @$n_ptr;

    my $ok = 1;
    my $span = 1;
    my $n = $nn[0];
    if($read < 3)
    {
	if($read == 1)
	{
	    $span = $nn[0];
	    $n = $nn[1]; 
	}
	elsif($read == 2)
	{
	    $span = $nn[0] * $nn[1];
	    $n = $nn[2]; 
	}

	$ok = 0;
	if(($nvals % $span == 0) && ($nvals < $n * $span))
	{
	    $ok = 1;
	}
    }
    if($ok == 1)
    {
	print OUT "$val\n";
    }
}
