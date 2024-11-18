#!/usr/bin/perl

my $narg = $#ARGV + 1;
if($narg != 2)
{
    print "usage: $0 <filename> <attribute name>\n";
    print "extracts points with values from vtk legacy file\n";
    exit();
}

my $filename = shift;
my $search_term = shift;
my $onpt = 0;
my $onply = 0;
my $onsclr = 0;
my $nel = 0;

my @ndx = ();
my @ndy = ();
my @ndz = ();
my @ptx = ();
my @pty = ();
my @ptz = ();

open(IN,"<$filename");
#print "open $filename\n";
while(<IN>)
{
    #section headers
    if(/POINTS\s(\S+)\s/)
    {
	#print "found POINTS header ... \n";
	$onpt = $1;
	$onply = 0;
	$onsclr = 0;
    }
    elsif(/POLYGONS\s(\S+)\s(\S+)/)
    {
	#print "found POLYGONS header ... \n";
	$onpt = 0;
	$onply = $1;
	$onsclr = 0;
	$nel = $1;
    }
    elsif(/SCALARS\s/)
    {
	#print "found SCALARS header ... \n";
	$onpt = 0;
	$onply = 0;
	$onsclr = 0;
	if($_=~m/$search_term/)
	{
	    my $line = <IN>;
	    $onsclr = 1;
	}
    }
    #evaluate flags
    if($onpt > 0)
    {
	#print "gettings $onpt pts ...\n";
	for(my $i=0;$i<$onpt;$i++)
	{
	    my $line = <IN>;
	    chomp($line);
	    if($line=~m/(\S+)\s(\S+)\s(\S+)/)
	    {
		push(@ndx,$1);
		push(@ndy,$2);
		push(@ndz,$3);
	    }
	}
	$onpt = 0;
    }
    elsif($onply > 0)
    {
	#print "gettings $onply plys ...\n";
	for(my $i=0;$i<$onply;$i++)
	{
	    my $line = <IN>;
	    chomp($line);
	    if($line=~m/\S+\s(\S+)\s(\S+)\s(\S+)/)
	    {
		#add to pts
		my $x = $ndx[$1];
		$x += $ndx[$2];
		$x += $ndx[$3];
		$x /= 3.0;
		my $y = $ndy[$1];
		$y += $ndy[$2];
		$y += $ndy[$3];
		$y /= 3.0;
		my $z = $ndz[$1];
		$z += $ndz[$2];
		$z += $ndz[$3];
		$z /= 3.0;
		#pts
		push(@ptx,$x);
		push(@pty,$y);
		push(@ptz,$z);
	    }
	}
	$onply = 0;
    }
    elsif($onsclr > 0)
    {
	#print "gettings $nel scalars ...\n";
	my $line;# = <IN>;
	for(my $i=0;$i<$nel;$i++)
	{
	    $line = <IN>;
	    chomp($line);
	    if($line=~m/^(\S+)/)
	    {
		print("$ptx[$i] $pty[$i] $ptz[$i] $1\n");
	    }
	}
	$onsclr = 0;
    }
}
close IN;
