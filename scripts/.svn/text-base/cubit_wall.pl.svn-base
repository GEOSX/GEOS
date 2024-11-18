#!/usr/bin/perl
use strict;
my $narg = $#ARGV + 1;
die "usage: $0 <minmax file> <element size> <wall thickness> <mesh output file name>\n" if($narg != 4);

my $filename = shift;
my $elsize = shift;
my $wallthickness = shift;
my $outfile = shift;

my @mm = ();

#read extrema
die "cannot open $filename\n" if(!open(IN,"<$filename"));
while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line=~s/^\s+//;
    my @tmp = split(/\s+/,$line);
    $tmp[0] = $tmp[2] - $tmp[1];
    push(@mm,\@tmp);
}
close IN;

#extract extrema and write
print "#!python\n";

my $ivol = 1;
my $cx = 0.5 * ($mm[0][1] + $mm[0][2]);
my $cy = 0.5 * ($mm[1][1] + $mm[1][2]);
my $cz = 0.5 * ($mm[2][1] + $mm[2][2]);
for(my $i = 0; $i < 3; $i++)
{
    my $dx = $mm[0][0]; 
    my $dy = $mm[1][0]; 
    my $dz = $mm[2][0];

    my $cdim = "";
    my $dd = 0;
    my $cc = 0;
    if($i == 0)
    {
	$cdim = "x";
	$dd = 0.5 * ($dx + $wallthickness);
	$cc = $cx;
	$dx = $wallthickness;
    }
    elsif($i == 1)
    {
	$cdim = "y";
	$dd = 0.5 * ($dy + $wallthickness);
	$cc = $cy;
	$dy = $wallthickness;
    }
    elsif($i == 2)
    {
	$cdim = "z";
	$dd = 0.5 * ($dz + $wallthickness);
	$cc = $cz;
	$dz = $wallthickness;
    }

    for(my $j = -1; $j < 2; $j += 2)
    {
	print "cubit.cmd(\'create brick x $dx y $dy z $dz\')\n";
	my $dds = $cc + $j * $dd;
	print "cubit.cmd(\'move volume $ivol $cdim $dds\')\n";
	$ivol++;
    }
}

#do other operations
print "cubit.cmd('volume 1 2 3 4 5 6 size $elsize')
cubit.cmd('volume 1 2 3 4 5 6 scheme Map')
cubit.cmd('mesh volume 1 2 3 4 5 6')
cubit.cmd('block 1 volume 1')
cubit.cmd('block 2 volume 2')
cubit.cmd('block 3 volume 3')
cubit.cmd('block 4 volume 4')
cubit.cmd('block 5 volume 5')
cubit.cmd('block 6 volume 6')
cubit.cmd('nodeset 1 volume 1')
cubit.cmd('nodeset 2 volume 2')
cubit.cmd('nodeset 3 volume 3');
cubit.cmd('nodeset 4 volume 4')
cubit.cmd('nodeset 5 volume 5')
cubit.cmd('nodeset 6 volume 6');
cubit.cmd('set Abaqus precision 6')
cubit.cmd('export Abaqus \"$outfile\" Block 1 2 3 4 5 6 Nodeset 1 2 3 4 5 6 dimension 3 everything overwrite cubitids ')
cubit.cmd('quit')
";
