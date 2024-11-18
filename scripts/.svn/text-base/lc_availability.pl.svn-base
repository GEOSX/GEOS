#!/usr/bin/perl
use strict;
#use File::HomeDir;
my $narg = $#ARGV + 1;
die "usage: $0 <min procs> <max procs> <min time> <max time> <weekday>\n" if($narg != 5);

#get arguments
my @procs = (0,0);
$procs[0] = shift;
$procs[1] = shift;
my @times = (0,0);
$times[0] = shift;
$times[1] = shift;
my $weekday = shift;

#get possible clusters
my %clusters = ();
#my $filename = File::HomeDir->my_home . "/_lcclusters";
my $filename = $ENV{"HOME"} . "/_lcclusters";
die "cannot open $filename" if(!open(IN,"<$filename"));
while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line=~s/^\s+//;
    my $ok = 0;
    if(/weekday/)
    {
	if($weekday==1)
	{
	    $ok = 1;
	}
    }
    else
    {
	if($weekday==0)
	{
	    $ok = 1;
	}
    }
    if($ok == 1)
    {
	my @arr = split(/\t/,$line);
	my $maxp = $arr[2]*$arr[3];
	my $maxt = $arr[4] * 3600;
	if($maxp >= $procs[0] && $maxt >= $times[0])
	{
	    my $cluster = $arr[0];
	    $clusters{$cluster} = \@arr;
	    print STDERR "-->showbf: $cluster eligible (i.e., has the required configuration and time constraints)\n";
	}
    }
}
close IN;
{
    my @keyarr = keys %clusters;
    if($#keyarr < 0)
    {
	exit(1);
    }
}

open(IN,"showbf|");
{
    while(<IN>)
    {
	my $line = $_;
	$line =~s/^\s+//;
	chomp($line);
	my @arr = split(/\s+/,$line);
	my $cluster = $arr[0];
	if(exists $clusters{$cluster})
	{
	    my $core_available = ${$clusters{$cluster}}[3] * $arr[2];
	    print STDERR "-->showbf: $cluster has $core_available cores now ($procs[0] requested)\n";
	    if($procs[0] <= $core_available)
	    {
		my $time_available = ${$clusters{$cluster}}[4] * 3600;
		print STDERR "-->showbf: ... but we need $time_available seconds\n";
		if($core_available > $procs[1])
		{
		    $core_available = $procs[1];
		}
		my $print_and_exit = 0;
		if(/INFINITY/)
		{
		    $print_and_exit = 1;
		    print STDERR "-->showbf: ... and we discovered it has as much time as we desire\n";
		}
		else
		{
		    my $time_available = 0;
		    my @tarr = split(/\:/,$arr[3]);
		    my $fct = 1;
		    $time_available += $tarr[$#tarr]*$fct;
		    $fct *= 60;
		    $time_available += $tarr[$#tarr-1]*$fct;
		    $fct *= 60;
		    $time_available += $tarr[$#tarr-2]*$fct;
		    $fct *= 24;
		    if($#tarr == 3)
		    {
			$time_available += $tarr[$#tarr-3]*$fct;
		    }
		    print STDERR "-->showbf: ... and we discovered it has $time_available seconds ($arr[3])\n";
		    if($time_available >= $times[0])
		    {
			$print_and_exit = 1;
		    }
		    else
		    {
			print STDERR "   so far no luck with showbf: $cluster is ready to go with $core_available cores and $time_available seconds\n";
		    }
		}
		if($print_and_exit == 1)
		{
		    my $nnodes = $core_available % ${$clusters{$cluster}}[3];
		    if($nnodes == 0)
		    {
			$nnodes = $core_available / ${$clusters{$cluster}}[3];
		    }
		    else
		    {
			$nnodes = ($core_available - $nnodes) / ${$clusters{$cluster}}[3] + 1;
		    }
		    print STDERR "We're in luck: $cluster is ready to go with $core_available cores and $time_available seconds\n";
		    print "$cluster\t$nnodes\t$core_available\t$time_available\n";
		    exit 0;
		}
	    }
	}
    }
}


#if there is not one immediately available, just grab cab or a random one
{
    my @keyarr = keys %clusters;
    my $cluster = $keyarr[0];
    my $time_available = ${$clusters{$cluster}}[4] * 3600;
    my $core_available = ${$clusters{$cluster}}[3] * ${$clusters{$cluster}}[2];
    if($core_available > $procs[1])
    {
	$core_available = $procs[1];
    }
    my $nnodes = $core_available % ${$clusters{$cluster}}[3];
    if($nnodes == 0)
    {
	$nnodes = $core_available / ${$clusters{$cluster}}[3];
    }
    else
    {
	$nnodes = ($core_available - $nnodes) / ${$clusters{$cluster}}[3] + 1;
    }
    print "$cluster\t$nnodes\t$core_available\t$time_available\n";
    exit 0;
}
