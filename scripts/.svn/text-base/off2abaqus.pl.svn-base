#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;

die "usage: $0 <filename(s)>\n" if($narg==0);
my @filenames = ();

for(my $i = 0; $i < $narg; $i++)
{
    my $filename = shift;
    push(@filenames, $filename);
}

my $ttt = gmtime();
print "*HEADING
Scott Johnson: $ttt
*NODE
";

my $ind = 1;
my @nnds = ();
my @els = ();
for(my $i = 0; $i <= $#filenames; $i++)
{
    my $filename = $filenames[$i];
    #die "$filename";
    open(IN,"<$filename");

    my $nel = 0;
    my $nnd = 0;
    my $offset = $ind;
    {
	my $l1 = <IN>;
	my $l2 = <IN>;
	chomp($l2);
	my @arr = split(/\s+/,$l2);
	$nnd = $arr[0];
	$nel = $arr[1];
    }
    push(@nnds,$nnd);

    for(my $i=0;$i<$nnd;$i++) {
	my $line = <IN>;
	chomp($line);
	$line=~s/^\s+//;
	my @arr = split(/\s+/,$line);
	#push(@nodes,[$ind, $arr[0],$arr[1],$arr[2]]);
	print "$ind,$arr[0],$arr[1],$arr[2]\n";
	$ind++;
    }

    for(my $i=0;$i<$nel;$i++)
    {
	my $line = <IN>;
	chomp($line);
	$line=~s/^\s+//;
	my @arr = split(/\s+/,$line);
	my @el = ();

#NOTE: QHULL PRODUCES CW NODE ORDERING
#      WE NEED CCW FOR GPAC, SO ...
	for(my $j = $#arr; $j > 0; $j--)
	{
	    push(@el, $arr[$j]+$offset);
	}
	push(@els,\@el);
    }

    close IN;
}

#print elements
print "*ELEMENT,TYPE=STRI3,ELSET=PM1\n";
my $iel = 1;
for(my $i = 0; $i <= $#els; $i++)
{
    my @loop = @{$els[$i]};
#    die "$#loop -- @loop\n" if($#loop != 2);
    for(my $j = 2; $j <= $#loop; $j++)
    {
	my $j0 = $j-1;
	my $j1 = $j;
	print "$iel,$loop[0],$loop[$j0],$loop[$j1]\n";
	$iel++;
    }
}

#print nodesets
$ind = 1;
for(my $i = 0; $i <= $#nnds; $i++)
{
    print "*NSET,NSET=de$i\n";
    for(my $j=0; $j<$nnds[$i]; $j++)
    {
	print "$ind";
	if($j==($nnds[$i]-1)) {
	    print "\n";
	} else {
	    print ",";
	}
	$ind++;
    }
}
