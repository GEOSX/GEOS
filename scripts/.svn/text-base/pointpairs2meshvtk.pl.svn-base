#!/bin/perl
use strict;
my $narg = $#ARGV + 1;
die "usage: <strike point pair file> <dip top-to-bottom depth-strike scale file> <center x, y, z> <scalar value at center> <scalar value at edge>\n" if($narg != 7);

my $fs = shift;
my $fd = shift;
my @x0 = (0,0,0);
for(my $i = 0; $i < 3; $i++)
{
	$x0[$i] = shift;
}
my $s0 = shift;
my $s1 = shift;

#LOAD STRIKE POINT PAIRS
my @cog = (0,0,0);	
my @ss = ();
{
	my $ddc = 0;
	die "cannot open $fs" if(!open(IN,"<$fs"));
	while(<IN>)
	{
		my @arr = &SplitMe($_);
		push(@ss, \@arr);
		my @pt0 = ($arr[0], $arr[1], 0);
		my @pt1 = ($arr[3], $arr[4], 0);
		my $dd = &Length(\@pt0, \@pt1);
		$ddc += $dd;
		$cog[0] += 0.5 * ($pt0[0] + $pt1[0]) * $dd;
		$cog[1] += 0.5 * ($pt0[1] + $pt1[1]) * $dd;
	}
	close IN;

	#get the cog and displacement
	$cog[0] /= $ddc;
	$cog[1] /= $ddc;
}

#get largest distance from cog in strike direction
my $dd_max = 0;
{
	my @pt = ($ss[0][0], $ss[0][1], 0);
	$dd_max = &Length(\@cog, $ss[0]);
	for(my $i = 0; $i <= $#ss; $i++)
	{
		$pt[0] = $ss[$i][3];
		$pt[1] = $ss[$i][4];
		my $dd = &Length(\@pt, \@cog);
		if($dd > $dd_max)
		{
			$dd_max = $dd;
		}
	}
}
my $ns = $#ss + 1;

#LOAD DIP-STRIKE_SCALE TUPLES
my @dss = ();
$cog[2] = 0;
{
	my $ddc = 0;
	die "cannot open $fd" if(!open(IN,"<$fd"));
	while(<IN>)
	{
		my @arr = &SplitMe($_);
		die "wrong number of entries\n" if($#arr != 1);
		if($#dss >= 0)
		{
			my $d0 = $dss[$#dss][0];
			my $d1 = $arr[0];
			my $ss0 = $dss[$#dss][1];
			my $ss1 = $arr[1];
			
			#get incremental depth
			my $ddepth = $d0 - $d1;
			die "depth must be specified in monotonically decreasing order with non-zero gradient!\n" if($ddepth <= 0);

			#increment area summation
			$ddc += $ddepth * 0.5 * ($ss0 + $ss1);

			#add COG components
			if($ss0 > $ss1)
			{
				#contribution of rectangular part
				$cog[2] += $ss1 * $ddepth * ($d1 + 0.5 * $ddepth);
				
				#contribution of triangular part
				$cog[2] += 0.5 * ($ss0 - $ss1) * $ddepth * ($d1 + (2/3) * $ddepth);
			}
			else
			{
				#contribution of rectangular part
				$cog[2] += $ss0 * $ddepth * ($d1 + 0.5 * $ddepth);
				
				#contribution of triangular part
				$cog[2] += 0.5 * ($ss1 - $ss0) * $ddepth * ($d1 + (1/3) * $ddepth);				
			}
		}
		
		#add current entry to the list
		push(@dss, \@arr);
	}
	$cog[2] /= $ddc;
	close IN;
}
my $nd = $#dss + 1;

#STORE NODES WITH SCALAR (AFTER RE-CENTERING)
my @nds = ();
for(my $j = 0; $j < $nd; $j++)
{
	my $d0 = $dss[$j][0];
	my $ss0 = $dss[$j][1];

	#get scalar value
	my $scalar = $s0;
	{
		my @tpt = ($ss[0][0], $ss[0][1], $d0);
		my $dd = &Length(\@tpt, \@cog);
		my $fct = $dd / $dd_max;
		if($fct < 1e-6)
		{
			$fct = 1e-6;
		}
		elsif($fct > 1.0)
		{
			$fct = 1.0;
		}
		$scalar += ($s1 - $s0) * $fct;
	}	
	my @pt0 = ($ss0 * ($ss[0][0] - $cog[0]) + $x0[0], $ss0 * ($ss[0][1] - $cog[1]) + $x0[1], $d0 + $x0[2] - $cog[2], $scalar);
	push(@nds,\@pt0);
	for(my $i = 0; $i < $ns; $i++)
	{
		$scalar = $s0;
		{
			my @tpt = ($ss[$i][3], $ss[$i][4], $d0);
			my $dd = &Length(\@tpt, \@cog);
			my $fct = $dd / $dd_max;
			if($fct < 1e-6)
			{
				$fct = 1e-6;
			}
			elsif($fct > 1.0)
			{
				$fct = 1.0;
			}
			$scalar += ($s1 - $s0) * $fct;
		}
		my @pt = ($ss0 * ($ss[$i][3] - $cog[0]) + $x0[0], $ss0 * ($ss[$i][4] - $cog[1]) + $x0[1], $d0 + $x0[2] - $cog[2], $scalar);
		push(@nds, \@pt);
	}
}

my $nnodes = ($ns + 1) * ($nd);
die "nnodes = $nnodes but length of node array is ".($#nds + 1)."\n" if($nnodes != ($#nds + 1));

my $npolys = ($ns * ($nd - 1));


#DEBUG
#{
#	for(my $j = 0; $j <= $#nds; $j++)
#	{
#		printf("%f %f %f\n", $nds[$j][0], $nds[$j][1], $nds[$j][2]);
#	}
#	die "debug";
#}

#print the file header
print "<?xml version=\"1.0\"?>\n";

#---------------------------------
#BEGIN TAGS
#---------------------------------
print "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n<PolyData>\n";
print "<Piece NumberOfPoints=\"$nnodes\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"$npolys\">\n";

#---------------------------------
#POINTS
#---------------------------------
print "<Points>\n<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
for(my $j = 0; $j <= $#nds; $j++)
{
	printf("%f %f %f\n", $nds[$j][0], $nds[$j][1], $nds[$j][2]);
}
print "</DataArray>\n</Points>\n";

#---------------------------------
#POLYS
#---------------------------------
print "<Polys>\n<DataArray type=\"UInt32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
my @poly_scalars = ();
for(my $j = 1; $j < $nd; $j++)
{
	for(my $i = 0; $i < $ns; $i++)
	{
		my $n00 = ($j - 1) * ($ns + 1) + $i;
		my $n10 = $n00 + 1;
		my $n01 = $n00 + ($ns + 1);
		my $n11 = $n01 + 1;
		print "$n00 $n10 $n11 $n01\n";
		my $scalar_interpolated = 0.25 * ($nds[$n00][3] + $nds[$n10][3] + $nds[$n01][3] + $nds[$n11][3]);
		push(@poly_scalars, $scalar_interpolated);
	}
}
print "</DataArray>\n";
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
{
    print "<DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"ascii\" Name=\"default\">\n";
	for(my $ii = 0; $ii <= $#poly_scalars; $ii++)
	{
		printf("%f\n", $poly_scalars[$ii]);
    }
    print "\n</DataArray>\n";
}
print "</CellData>\n";

#---------------------------------
#END TAGS
#---------------------------------
print "</Piece>\n</PolyData>\n</VTKFile>\n";

sub SplitMe
{
	my $line = shift;
	chomp($line);
	$line =~s/^\s+//;
	$line =~s/\,/ /;
	$line =~s/\s+/ /;
	my @arr = split(/\s+/,$line);
	return @arr;
}

sub Length
{
	my $ptr0 = shift;
	my $ptr1 = shift;
	my @arr0 = @$ptr0;
	my @arr1 = @$ptr1;
	my $dd = 0;
	for(my $i = 0; $i < 3; $i++)
	{
		my $d = $arr1[$i] - $arr0[$i];
		$d *= $d;
		$dd += $d;
	}
	$dd = sqrt($dd);
	return $dd;
}

sub InterpolatePoint
{
	my $ptr0 = shift;
	my $ptr1 = shift;
	my @arr0 = @$ptr0;
	my @arr1 = @$ptr1;
	my $fct = shift;
	
	my @ret = ();
	for(my $i = 0; $i < 3; $i++)
	{
		$ret[$i] = $arr0[$i] + ($arr1[$i] - $arr0[$i]) * $fct;
	}
	return @ret;
}



