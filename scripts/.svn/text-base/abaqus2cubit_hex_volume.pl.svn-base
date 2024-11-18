#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;

die "usage: $0 <element size> <Abaqus filename in> <Abaqus filename out>\nNote: this prints an appropriate Python file for CUBIT to STDOUT\n" if($narg!=3);

#read input
my $elsize = shift;
my $ifile = shift;
my $ofile = shift;

#print resultant CUBIT file

print "#!python
cubit.cmd('import abaqus mesh geometry  \"$ifile\"')
cubit.cmd('create volume 1 surface 1 to 1000')
cubit.cmd('volume 1 to 1000 size $elsize')
cubit.cmd('volume 1 to 1000 scheme Map')
cubit.cmd('mesh volume 1 to 1000')
cubit.cmd('set duplicate block elements off')
cubit.cmd('block 1 volume 1 to 1000 ')
cubit.cmd('set Abaqus precision 6')
cubit.cmd('export Abaqus \"$ofile\" Block 1 dimension 3 overwrite')
cubit.cmd('exit')
";
