#!/usr/bin/perl
use strict;
use Math::Trig;

my $narg = $#ARGV + 1;
my $toggle = 0;
if($toggle)
{
    die "usage: $0 <ABAQUS input filename> <fault definition file: space-delimited tuples of [nx,ny,nz,ux,uy,uz,lu,lv,x,y,z]>\n" if($narg != 2);
}
else
{
    die "usage: $0 <ABAQUS input filename> <fault definition file: space-delimited tuples of [x,y,z,length,depth,dip(deg),strike(deg)]>\n" if($narg != 2);
}

my $file = shift;
my $ffile = shift;

#---------------------------
#  GET FAULTS
#---------------------------
my @faults = ();
die "unable to open $ffile" if(!open(IN,"<$ffile"));
open(LOG,">log");
if($toggle)
{
    while(<IN>)
    {
	my $line = $_;
	$line =~s/^\s+//;
	chomp($line);
	my @arr = split(/\s+/,$line);
	
	die "fault entry must be length 11: nx ny nz ux uy uz lu lv x y z!!" if($#arr != 10);
	
	my $dot = $arr[0]*$arr[3] + $arr[1]*$arr[4] + $arr[2]*$arr[5];
	die "normal and u must be orthogonal!!" if($dot < -1e-6 || $dot > 1e-6);
	
	my $nn = $arr[0]*$arr[0] + $arr[1]*$arr[1] + $arr[2]*$arr[2];
	my $uu = $arr[3]*$arr[3] + $arr[4]*$arr[4] + $arr[5]*$arr[5];
	die "normal must be unit vector!!" if(($nn-1) < -1e-6 || ($nn-1) > 1e-6);
	die "in plane vector must be unit vector!!" if(($uu-1) < -1e-6 || ($uu-1) > 1e-6);
	
	push(@faults,\@arr);
    }
}
else
{
    my $pi = 4 * atan(1);
    open(OUT,">pp_faults");
    while(<IN>)
    {
	my $line = $_;
	$line =~s/^\s+//;
	chomp($line);
	my @arr = split(/\s+/,$line);
	
	die "fault entry must be length 7: x y z length depth dip strike!!" if($#arr != 6);
	
	#get values
	my $x = $arr[0];
	my $y = $arr[1];
	my $z = $arr[2];
	my $ll = $arr[3];
	my $depth = $arr[4];
	my $dip = $arr[5];
	my $strike = $arr[6];

	#calculate normal
	my $nz = cos($dip*$pi/180.0);
	my $nxy = sqrt(1.0 - $nz * $nz);
	my $nx = $nxy * cos($strike*$pi/180.0);
	my $ny = -$nxy * sin($strike*$pi/180.0);
	
	#get ux uy uz (vector along strike)
	my $ux = sin($strike*$pi/180.0);
	my $uy = cos($strike*$pi/180.0);
	my $uz = 0;

	my $x0 = $x - $ux * $ll * 0.5;
	my $x1 = $x + $ux * $ll * 0.5;
	my $y0 = $y - $uy * $ll * 0.5;
	my $y1 = $y + $uy * $ll * 0.5;

	print OUT "$x0 $y0 0 $x1 $y1 0\n";

	#print LOG "$nx $ny $nz $ux $uy $uz $ll $depth $x $y $z\n";

	my @arr2 = ($nx,$ny,$nz,$ux,$uy,$uz,$ll,$depth,$x,$y,$z);
	push(@faults,\@arr2);
    }
    close IN;
    close OUT;
}
my $nfaults = $#faults;
{
    $nfaults++; 
    print LOG "READ $nfaults FAULTS\n";
} 
close IN;

#---------------------------
#  PRINT E.T. AND GET NODES/ELEMENTS
#---------------------------
die "unable to open $file" if(!open(IN,"<$file"));

my $node_on = 0;
my $element_on = 0;
my $endstr = "";
my $end_on = 0;

my @nodes = ();
my @elements = ();
my @faces = ([0,3,2,1],[4,5,6,7],[0,4,7,3],[1,2,6,5],[2,3,7,6],[0,1,5,4]);
while(<IN>)
{
    if($end_on == 1)
    {
	$endstr .= $_;
    }
    elsif(/^\*(\S+)/)
    {
	print $_;
	if(/NODE/)
	{
	    $node_on = 1;
	}
	else
	{
	    $node_on=0;
	}
	if(/PROPERTIES/)
	{
	    $endstr .= $_;
	    $end_on = 1;
	}
	elsif(/ELEMENT/)
	{
	    $element_on = 1;
	}
	else
	{
	    $element_on = 0;
	}
    }
    else
    {
	print $_;
	my $line = $_;
	chomp($line);
	my @arr = split(/\,\s+/,$line);
	
	if($node_on)
	{
	    $nodes[$arr[0]] = [$arr[1],$arr[2],$arr[3]];
	}
	elsif($element_on)
	{
	    push(@elements,[$arr[1],$arr[2],$arr[3],$arr[4],$arr[5],$arr[6],$arr[7],$arr[8]]);
	}
    }
}
my $nnds = $#nodes + 1;
print LOG "READ $nnds NODES\n";
my $nels = $#elements + 1;
print LOG "READ $nels ELEMENTS\n";
close IN;

#---------------------------
#  INITIALIZE THE NEW NODESETS
#---------------------------

my @nf = ();
for(my $i = 0; $i <= $#nodes; $i++)
{
    $nf[$i] = 0;
}

#---------------------------
#  CALCULATE THE NEW NODESETS
#---------------------------

#AN ELEMENT MAY CONTRIBUTE NODES IF AT LEAST ONE OF ITS NODES IS ON EACH SIDE OF THE PLANE.
#OF THESE NODES, IF A FACE HAS ALL OF ITS NODES AS BEING FLAGGED ELIGIBLE THEN ALL OF ITS
#NODES ARE ADDED TO THE NODESET.
for(my $i = 0; $i <= $#elements; $i++)
{
    for(my $ifault = 0; $ifault < $nfaults; $ifault++)
    {
	my @fnrm = ($faults[$ifault][0],$faults[$ifault][1],$faults[$ifault][2]);
	my @fu = ($faults[$ifault][3],$faults[$ifault][4],$faults[$ifault][5]);
	my @fv = ($fnrm[1]*$fu[2] - $fnrm[2]*$fu[1], $fnrm[2]*$fu[0] - $fnrm[0]*$fu[2], $fnrm[0]*$fu[1] - $fnrm[1]*$fu[0]);
	my $lu = 0.5 * $faults[$ifault][6];
	my $lv = 0.5 * $faults[$ifault][7];
	my @fx = ($faults[$ifault][8],$faults[$ifault][9],$faults[$ifault][10]);
	my $npos = 0;
	my $nneg = 0;
	my @inodes = ();
	for(my $j = 0; $j < 8; $j++)
	{
	    my $inode = $elements[$i][$j];
	    my $dd = 0;
	    my $uu = 0;
	    my $vv = 0;
	    for(my $k = 0; $k < 3; $k++)
	    {
		my $xx = $nodes[$inode][$k]-$fx[$k];
		$dd += $xx * $fnrm[$k];
		$uu += $xx * $fu[$k];
		$vv += $xx * $fv[$k];
	    }
	    if($dd >= 0)
	    {
		$npos = 1;
		if($uu <= $lu && $uu >= -$lu && $vv <= $lv && $vv >= -$lv)
		{
		    push(@inodes,$inode);
		}
	    }
	    else
	    {
		$nneg = 1;
	    }
	}
	#for each node on the element

	#if the nodes are on opposite sides of the plane
	#take all those on the positive side and put in the 
	#prospective node list
	if(($nneg * $npos) == 1)
	{
	    for(my $j = 0; $j <= $#inodes; $j++)
	    {
		$nf[$inodes[$j]] = 1;
	    }
	}
    }
    #for each fault
}
#for each element

open(OUT,">pp");
for(my $i = 0; $i <= $#elements; $i++)
{
    #now, go through all faces on the element
    for(my $j = 0; $j <= $#faces; $j++)
    {
	my $prod = 1;
	for(my $k = 0; $k < 4; $k++)
	{
	    my $inode = $elements[$i][$faces[$j][$k]];
	    $prod *= $nf[$inode];
	}
	if($prod > 0)
	{
	    for(my $k = 0; $k < 4; $k++)
	    {
		my $inode0 = $elements[$i][$faces[$j][$k]];
		my $inode1 = 0;
		if($k == 3)
		{
		    $inode1 = $elements[$i][$faces[$j][0]];
		}
		else
		{
		    $inode1 = $elements[$i][$faces[$j][$k+1]];
		}
		print OUT "$nodes[$inode0][0] $nodes[$inode0][1] $nodes[$inode0][2] $nodes[$inode1][0] $nodes[$inode1][1] $nodes[$inode1][2]\n";
		$nf[$inode0] = 2;
	    }
	} 
    }
}
close OUT;

#---------------------------
#  PRINT THE NEW NODESETS
#---------------------------
print "*NSET, NSET=FF0\n"; 
my $ii = 0;
my $ipts = 0;
open(OUT,">npts");
for(my $i = 0; $i <= $#nf; $i++)
{
    if($nf[$i]==2)
    {
	print OUT "$nodes[$i][0] $nodes[$i][1] $nodes[$i][2]\n";
	my $str = sprintf("%06d,",$i);
	if($ii==6)
	{
	    print "$str\n";
	    $ii = 0;
	}
	else
	{
	    print "$str";
	    $ii++;
	}
	$ipts++;
    }
}
close OUT;
print "**\n";

print LOG "FOUND $ipts NODES IN THE FLOW FACES COLLECTION\n";
#---------------------------
#  PRINT THE LAST STUFF
#---------------------------
print $endstr;
close LOG;
exit();





















#	for(my $j = 0; $j <= $#faces; $j++)
#	{
#	    my @n = (0,0,0);
#	    #cross 1
#	    {
#		my @u = (0,0,0);
#		my @v = (0,0,0);
#		for(my $jj = 0; $jj < 3; $jj++)
#		{
#		    $u[$jj] = $nodes[$elements[$i][$faces[$j][1]]][$jj] - $nodes[$elements[$i][$faces[$j][0]]][$jj];
#		    $v[$jj] = $nodes[$elements[$i][$faces[$j][2]]][$jj] - $nodes[$elements[$i][$faces[$j][0]]][$jj];
#		}
#		$n[0] += $u[1]*$v[2] - $u[2]*$v[1];
#		$n[1] += $u[2]*$v[0] - $u[0]*$v[2];
#		$n[2] += $u[0]*$v[1] - $u[1]*$v[0];
#	    }
#	    #cross 2
#	    {
#		my @u = (0,0,0);
#		my @v = (0,0,0);
#		for(my $jj = 0; $jj < 3; $jj++)
#		{
#		    $u[$jj] = $nodes[$elements[$i][$faces[$j][2]]][$jj] - $nodes[$elements[$i][$faces[$j][0]]][$jj];
#		    $v[$jj] = $nodes[$elements[$i][$faces[$j][3]]][$jj] - $nodes[$elements[$i][$faces[$j][0]]][$jj];
#		}
#		$n[0] += $u[1]*$v[2] - $u[2]*$v[1];
#		$n[1] += $u[2]*$v[0] - $u[0]*$v[2];
#		$n[2] += $u[0]*$v[1] - $u[1]*$v[0];
#	    }
#	    #normalize
#	    {
#		my $sum = 0;
#		for(my $jj = 0; $jj < 3; $jj++)
#		{
#		    $sum += $n[$jj]*$n[$jj];
#		}
#		if($sum > 0)
#		{
#		    $sum = 1.0 / sqrt($sum);
#		    for(my $jj = 0; $jj < 3; $jj++)
#		    {
#			$n[$jj] *= $sum;
#		    }
#		}
#	    }
	    
	    
	    
