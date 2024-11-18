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
    print "usage: $0 <input_file> <scalar0> <scalar1> ... <scalarN>\n";
    print "input is space-delimited, carriage return terminated\n";
    print "note: positions must be defined by the field names 'x', 'y', and 'z'\n";
    exit();
}

my @args = ();
my $nscalars = $narg-1;
my @pos = (-1,-1,-1);
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
	if($arg eq "x"){$pos[0]=$i;}
	if($arg eq "y"){$pos[1]=$i;}
	if($arg eq "z"){$pos[2]=$i;}
	push(@args,$arg);
    }
}

#print the file header
print "<?xml version=\"1.0\"?>\n";
print "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n";

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

print "<Piece NumberOfPoints=\"$npoint\" NumberOfCells=\"$npoint\">\n";
print "<Points>\n<DataArray type=\"Float32\" Name=\"SPH_Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";

my $ndim = 2;
if($pos[2]>-1)
{
    $ndim = 3;
}

#print "$pos[0] $pos[1] $pos[2] $ndim\n";exit();

#print points
{
    open(IN,"<$file");
    while(<IN>)
    {
	#get position
	my $line = $_;chomp($line);
	my @values = split(/\s+/,$line);
	for(my $i=0;$i<$ndim;$i++)
	{
	    my $value = $values[$pos[$i]];
	    $value = unpack("f",pack("f",$value));
	    print sprintf('%f',$value);
	    print " ";
	}
	for(my $i=$ndim;$i<3;$i++)
	{
	    print "0 ";
	}
	#exit();
    }
    close IN;
}
print "\n</DataArray>\n</Points>\n<Cells>\n";

#need connectivity for some reason for this to work; start with 0-based for "connectivity"
print "<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\" NumberOfComponents=\"1\">\n";
for(my $i=0;$i<$npoint;$i++)
{
    print "$i ";
}
print "\n</DataArray>\n";

#then with 1-based for "offsets"
print "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\" NumberOfComponents=\"1\">\n";
for(my $i=1;$i<=$npoint;$i++)
{
    print "$i ";
}
print "\n</DataArray>\n";

#then with 1's for "types"
print "<DataArray type=\"UInt32\" Name=\"types\" format=\"ascii\" NumberOfComponents=\"1\">\n";
for(my $i=1;$i<=$npoint;$i++)
{
    print "1 ";
}
print "\n</DataArray>\n";
print "</Cells>\n<PointData>\n";

#start scalars offset by 3 (position)
for(my $i=0;$i<$nscalars;$i++)
{
    if($i!=$pos[0] && $i!=$pos[1] && $i!=$pos[2])
    {
	my $sname = $args[$i];

	#print each scalar and the associated data ... make sure it is Float32!!
	print "<DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"ascii\" Name=\"$sname\">\n";
	#print scalars
	{
	    open(IN,"<$file");
	    while(<IN>)
	    {
		#get position
		my $line = $_;chomp($line);
		my @values = split(/\s+/,$line);
		my $value = $values[$i];
		#print "$value $i\n";
		#exit();
		$value = unpack("f",pack("f",$value));
		print sprintf('%f',$value);
		print " ";
	    }
	    close IN;
	}
	print "\n</DataArray>\n";
    }
}
print "</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
