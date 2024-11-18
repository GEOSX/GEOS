#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <filename>\n" if($narg != 1);

my $file = shift;
open(IN,"<$file");
my $node_on = 0;
my $el_on = 0;

my @els = ();

my $now = localtime time;

print "*HEADING
Tied nodes generated from $file at $now
**
********************************** P A R T S **********************************
*PART, NAME=Part-Default
**
********************************** N O D E S **********************************
*NODE, NSET=ALLNODES
";
{
    my $nodenum = 1;
    my $node_offset = 0;
    my $elnum = 1;
    my @nodes = ();
    my @mask = (0,4,6,2,1,5,7,3);
    while(<IN>)
    {
	if(/^\*(\S+)/)
	{
	    if(/NODE/)
	    {
		$node_on = 1;
	    }
	    else
	    {
		$node_on=0;
	    }
	    if(/ELEMENT/)
	    {
		if(/TYPE/)
		{
		    $el_on = 1;
		}
		else
		{
		    $el_on=0;
		}
	    }
	    else
	    {
		$el_on=0;
	    }
	}
	else
	{
	    if($node_on == 1 || $el_on == 1)
	    {
		my $line = $_;
		chomp($line);
		$line=~s/^s+//;
		my @arr = split(/\,\s*/,$line);

		if($node_on == 1)
		{
		    my @xnd = ($arr[1],$arr[2],$arr[3]);
		    $nodes[$arr[0]] = \@xnd;
		    #die "$nodes[$arr[0]][0] $nodes[$arr[0]][1] $nodes[$arr[0]][2]\n";
		}
		elsif($el_on == 1)
		{
		    my @local_nodes = ();

		    #die "$arr[0]: $arr[1] $arr[2] $arr[3] $arr[4] $arr[5] $arr[6] $arr[7] $arr[8]: $line";
		    
                    #get local nodes
		    for(my $i = 1; $i <= $#arr; $i++)
		    {
			my @xnd = ($nodes[$arr[$i]][0], $nodes[$arr[$i]][1], $nodes[$arr[$i]][2]);
			#die "$xnd[0] $xnd[1], $xnd[2]\n" if($i==2);
			$local_nodes[$mask[$i-1]] = \@xnd;
		    }

		    #print nodes in order
		    for(my $i = 1; $i <= 27; $i++)
		    {
			print "$nodenum, ";
			&PrintNodalLocation ($i, \@local_nodes);
			$nodenum++;
		    }

		    #create the octants
		    
		    # NUMBER 1
		    {
			my @el = ($elnum,1,10,13,4,2,11,14,5);
			for(my $i = 1; $i <= $#el; $i++)
			{
			    $el[$i] += $node_offset;
			}
			push(@els,\@el);
			$elnum++;
		    }

		    # NUMBER 2
		    {
			my @el = ($elnum,2,11,14,5,3,12,15,6);
			for(my $i = 1; $i <= $#el; $i++)
			{
			    $el[$i] += $node_offset;
			}
			push(@els,\@el);
			$elnum++;
		    }

		    # NUMBER 3
		    {
			my @el = ($elnum,4,13,16,7,5,14,17,8);
			for(my $i = 1; $i <= $#el; $i++)
			{
			    $el[$i] += $node_offset;
			}
			push(@els,\@el);
			$elnum++;
		    }

		    # NUMBER 4
		    {
			my @el = ($elnum,5,14,17,8,6,15,18,9);
			for(my $i = 1; $i <= $#el; $i++)
			{
			    $el[$i] += $node_offset;
			}
			push(@els,\@el);
			$elnum++;
		    }


		    # NUMBER 5
		    {
			my @el = ($elnum,10,19,22,13,11,20,23,14);
			for(my $i = 1; $i <= $#el; $i++)
			{
			    $el[$i] += $node_offset;
			}
			push(@els,\@el);
			$elnum++;
		    }

		    # NUMBER 6
		    {
			my @el = ($elnum,13,22,25,16,14,23,26,17);
			for(my $i = 1; $i <= $#el; $i++)
			{
			    $el[$i] += $node_offset;
			}
			push(@els,\@el);
			$elnum++;
		    }

		    # NUMBER 7
		    {
			my @el = ($elnum,11,20,23,14,12,21,24,15);
			for(my $i = 1; $i <= $#el; $i++)
			{
			    $el[$i] += $node_offset;
			}
			push(@els,\@el);
			$elnum++;
		    }

		    # NUMBER 8
		    {
			my @el = ($elnum,14,23,26,17,15,24,27,18);
			for(my $i = 1; $i <= $#el; $i++)
			{
			    $el[$i] += $node_offset;
			}
			push(@els,\@el);
			$elnum++;
		    }

		    $node_offset += 27;
		}
	    }
	}
    }
}
close IN;

print "**** WARNING - MATERIAL    1 UNDEFINED FOR THE FOLLOWING ELEMENTS ****
*ORIENTATION,NAME=SOR1,DEFINITION=COORDINATES
0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00
0,0.000000E+00
*SOLID SECTION,ELSET=PM1,MATERIAL=unknown,ORIENTATION=SOR1
*ELEMENT,TYPE=C3D8,ELSET=PM1
";

for(my $i = 0; $i <= $#els; $i++)
{
    my @arr = @{$els[$i]};
    for(my $j = 0; $j <= $#arr; $j++)
    {
	print "$els[$i][$j]";
	if($j < $#arr)
	{
	    print ",";
	}
	else
	{
	    print "\n";
	}
    }
}

#print "*NSET,NSET=block1
#1,5,3,7,2,6,4,8
#*NSET,NSET=block2
#9,13,11,15,10,14,12,16
#*NSET,NSET=xneg
#9,10,11,12
#*NSET,NSET=xpos
#5,6,7,8
#";


sub PrintNodalLocation
{
    my $local_node_number = shift;
    my $local_nodes_pt = shift;
    my @local_nodes = @$local_nodes_pt;

    die "out of range" if($local_node_number < 1 || $local_node_number > 27);

    my @x = (0,0,0);
    if($local_node_number == 1)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = $local_nodes[0][$i];
	}
    }
    elsif($local_node_number == 2)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[0][$i] + $local_nodes[1][$i]);
	}
    }
    elsif($local_node_number == 3)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = $local_nodes[1][$i]
	}
    }
    elsif($local_node_number == 4)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[0][$i] + $local_nodes[2][$i]);
	}
    }
    elsif($local_node_number == 5)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[0][$i] + $local_nodes[3][$i]);
	}
    }
    elsif($local_node_number == 6)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[1][$i] + $local_nodes[3][$i]);
	}
    }
    elsif($local_node_number == 7)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = $local_nodes[2][$i];
	}
    }
    elsif($local_node_number == 8)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[2][$i] + $local_nodes[3][$i]);
	}
    }
    elsif($local_node_number == 9)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = $local_nodes[3][$i];
	}
    }


    elsif($local_node_number == 10)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[0][$i] + $local_nodes[4][$i]);
	}
    }
    elsif($local_node_number == 11)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[0][$i] + $local_nodes[5][$i]);
	}
    }
    elsif($local_node_number == 12)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[1][$i] + $local_nodes[5][$i]);
	}
    }
    elsif($local_node_number == 13)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[0][$i] + $local_nodes[6][$i]);
	}
    }
    elsif($local_node_number == 14)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[3][$i] + $local_nodes[4][$i]);
	}
    }
    elsif($local_node_number == 15)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[3][$i] + $local_nodes[5][$i]);
	}
    }
    elsif($local_node_number == 16)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[2][$i] + $local_nodes[6][$i]);
	}
    }
    elsif($local_node_number == 17)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[2][$i] + $local_nodes[7][$i]);
	}
    }
    elsif($local_node_number == 18)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[3][$i] + $local_nodes[7][$i]);
	}
    }


    elsif($local_node_number == 19)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = $local_nodes[4][$i];
	}
    }
    elsif($local_node_number == 20)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[4][$i] + $local_nodes[5][$i]);
	}
    }
    elsif($local_node_number == 21)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = $local_nodes[5][$i];
	}
    }
    elsif($local_node_number == 22)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[4][$i] + $local_nodes[6][$i]);
	}
    }
    elsif($local_node_number == 23)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[4][$i] + $local_nodes[7][$i]);
	}
    }
    elsif($local_node_number == 24)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[5][$i] + $local_nodes[7][$i]);
	}
    }
    elsif($local_node_number == 25)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = $local_nodes[6][$i];
	}
    }
    elsif($local_node_number == 26)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = 0.5 * ($local_nodes[6][$i] + $local_nodes[7][$i]);
	}
    }
    elsif($local_node_number == 27)
    {
	for(my $i = 0; $i < 3; $i++)
	{
	    $x[$i] = $local_nodes[7][$i];
	}
    }


    #die "($local_nodes[0][0] $local_nodes[0][1] $local_nodes[0][2]) + ($local_nodes[1][0] $local_nodes[1][1] $local_nodes[1][2])\n";

    print "$x[0], $x[1], $x[2]\n";    
}
