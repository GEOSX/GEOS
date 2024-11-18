#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;

#determine if the user is asking for help (either "help" or no arguments)
my $print_help = 0;
my $file = "";
if($narg > 0)
{
    $file = shift;
    if($file eq "help") {$print_help=1;}
}
else
{
    $print_help = 1;
}

if($print_help)
{
    print "usage: $0 <file>\n";
    print "input is space-delimited, carriage return terminated\n";
    print "note: file format is <x,y,z> <nx,ny,nz> <ux,uy,uz> <du> <dv> <scalar 0 .. N>\n";
    exit();
}

my $nscalars = 0;

#get number of scalars
{
    open(IN,"<$file");
    my $line = <IN>;
    chomp($line);
    $line=~s/^\s+//;
    close IN;
    my @arr = split(/\s+/,$line);
    $nscalars = $#arr - 10;
}

#go through the file and render
my @planes = ();
{
    open(IN,"<$file");
    while(<IN>)
    {
	my $line = $_;
	chomp($line);
	$line=~s/^\s+//;
	my @arr = split(/\s+/,$line);
	push(@planes, \@arr);
    }
    close IN;
}

#print the file header
print "<?xml version=\"1.0\"?>\n";

#---------------------------------
#BEGIN TAGS
#---------------------------------
print "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n<PolyData>\n";

#count points
my $nnodes = 4*($#planes + 1);
my $npolys = $#planes + 1;
print "<Piece NumberOfPoints=\"$nnodes\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"$npolys\">\n";

#---------------------------------
#POINTS
#---------------------------------
print "<Points>\n<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
{
    open(IN,"<$file");
    my $inx = 3;
    my $iux = 6;
    my $ix = 0;
    my $idu = 9;
    my $idv = 10;
    while(<IN>)
    {
	my $line = $_;
	chomp($line);
	$line=~s/^\s+//;
	my @arr = split(/\s+/,$line);

	my $du = $arr[$idu];
	my $dv = $arr[$idv];
	my @xx = ($arr[$ix],$arr[$ix+1],$arr[$ix+2]);
	my @n = ($arr[$inx],$arr[$inx+1],$arr[$inx+2]);
	my @u = ($arr[$iux],$arr[$iux+1],$arr[$iux+2]);
	my @v = (0,0,0);

	#die "$u[0] $u[1] $u[2] $v[0] $v[1] $v[2]\n";

	#n cross u = v
	$v[0] = $n[1]*$u[2] - $u[1]*$n[2];
	$v[1] = $n[2]*$u[0] - $u[2]*$n[0];
	$v[2] = $n[0]*$u[1] - $u[0]*$n[1];

	my @nodes = ();
	#-u,-v; -u,v; u,-v; u,v
	for(my $i = -1; $i < 2; $i+=2)
	{
	    for(my $j = -1; $j < 2; $j+=2)
	    {
		for(my $ii = 0; $ii < 3; $ii++)
		{
		    my $x = $xx[$ii] + $i*$du*$u[$ii] + $j*$dv*$v[$ii]; 
		    print sprintf("%f ",$x);
		}
	    }
	}
    }
    close IN;
}
print "\n</DataArray>\n</Points>\n";

#---------------------------------
#POLYS
#---------------------------------
print "<Polys>\n<DataArray type=\"UInt32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
my @cxn = (0,2,3,1);
my $ie = 0;
for(my $i = 0; $i < $npolys; ++$i)
{
    for(my $ii = 0; $ii <= $#cxn; ++$ii)
    {
	my $value = $cxn[$ii] + $ie;
	print "$value ";
    }
    $ie += $#cxn + 1;
}
print "\n</DataArray>\n";
print "<DataArray type=\"UInt32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
for(my $i = 1; $i <= $npolys; $i++)
{
    my $value = 4 * $i;
    print "$value ";
}
print "\n</DataArray>\n</Polys>\n";

#---------------------------------
#CELLS
#---------------------------------
print "<CellData>\n";
#start scalars offset by 3 (position)
for(my $i=0;$i<$nscalars;$i++)
{
    my $sname = "s$i";
    my $iii = 10 + $i;
    
    #print each scalar and the associated data ... make sure it is Float32!!
    print "<DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"ascii\" Name=\"$sname\">\n";
    #print scalars
    open(IN,"<$file");
    while(<IN>)
    {
	my $line = $_;chomp($line);
	$line =~s/^\s+//;
	my @values = split(/\s+/,$line);
	my $value = $values[$iii];
	#print "$value $i\n";
	#exit();
	$value = unpack("f",pack("f",$value));
	print sprintf('%7.5e',$value);
	print " ";
    }
    close IN;
    print "\n</DataArray>\n";
}
print "</CellData>\n";

#---------------------------------
#END TAGS
#---------------------------------
print "</Piece>\n</PolyData>\n</VTKFile>\n";
