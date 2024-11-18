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
    print "usage: $0 <input_file> <dx> <scalar0> <scalar1> ... <scalarN>\n";
    print "input is space-delimited, carriage return terminated\n";
    print "note: positions must be defined by the field names 'x', 'y', and 'z'\n";
    exit();
}

my $dx = shift;
my @args = ();
my $nscalars = $narg-2;
my @pos = (-1,-1,-1);
my @sindices = ();
if($nscalars==0)
{
    #by default, assume the file only contains points
    push(@args,"x");
    push(@args,"y");
    push(@args,"z");
}
else
{
    for(my $i=0;$i<$nscalars;$i++)
    {
	my $arg = shift;
	if($arg eq "x"){$pos[0]=$i;next;}
	if($arg eq "y"){$pos[1]=$i;next;}
	if($arg eq "z"){$pos[2]=$i;next;}
	push(@args,$arg);
	push(@sindices,$i);
    }
}
$nscalars = $#args + 1;

#print the file header
print "<?xml version=\"1.0\"?>\n";

#---------------------------------
#BEGIN TAGS
#---------------------------------
print "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n<PolyData>\n";

#count points
my $npoint = 0;
{
    open(IN,"<$file");
    while(<IN>)
    {
	$npoint++;
    }
    close IN;
}
#print "npt=$npoint\n";
#exit();

my $nnodes = 8 * $npoint;
my $npolys = 6 * $npoint;
print "<Piece NumberOfPoints=\"$nnodes\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"$npolys\">\n";

#---------------------------------
#POINTS
#---------------------------------
print "<Points>\n<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
my $ndim = 2;
if($pos[2]>-1)
{
    $ndim = 3;
}
#print points
{
    my @nodes = (1,1,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,0,1,1,1,1,1);
    open(IN,"<$file");
    while(<IN>)
    {
	#get position
	my @pt = (0,0,0);
	{
	    my $line = $_;chomp($line);
	    my @values = split(/\s+/,$line);
	    for(my $i=0;$i<$ndim;$i++)
	    {
		my $value = $values[$pos[$i]];
		$value = unpack("f",pack("f",$value));
		$pt[$i] = $value;
	    }
	    for(my $i=$ndim;$i<3;$i++)
	    {
		$pt[$i] = 0;
	    }
	}
	#print positions of nodes
	for(my $i = 0; $i < 8; $i++)
	{
	    for(my $ii = 0; $ii < 3; $ii++)
	    {
		my $value = ($nodes[$i*3 + $ii] - 1) * $dx + $pt[$ii];
		print sprintf('%f',$value);
		print " ";
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
my @cxn = (0,1,2,3,2,1,4,5,3,2,5,6,1,0,7,4,0,3,6,7,5,4,7,6);
for(my $i = 0; $i < $npoint; $i++)
{
    for(my $ii = 0; $ii <= $#cxn; $ii++)
    {
	my $value = $cxn[$ii] + 8*$i;
	print "$value ";
    }
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
    my $sname = $args[$i];
    my $iii = $sindices[$i];
    
    #print each scalar and the associated data ... make sure it is Float32!!
    print "<DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"ascii\" Name=\"$sname\">\n";
    #print scalars
    open(IN,"<$file");
    while(<IN>)
    {
	#get position
	my $line = $_;chomp($line);
	my @values = split(/\s+/,$line);
	my $value = $values[$iii];
	#print "$value $i\n";
	#exit();
	$value = unpack("f",pack("f",$value));
	for(my $ii = 0; $ii < 6; $ii++)
	{
	    print sprintf('%7.5e',$value);
	    print " ";
	}
    }
    close IN;
    print "\n</DataArray>\n";
}
print "</CellData>\n";

#---------------------------------
#END TAGS
#---------------------------------
print "</Piece>\n</PolyData>\n</VTKFile>\n";
