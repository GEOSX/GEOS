#!/usr/bin/perl
my $narg = $#ARGV + 1;
die "usage: $0 <filename> <factor>\n" if($narg != 2);

my $filename = shift;
my $factor = shift;

die "cannot open $filename\n" if(!open(IN,"<$filename"));
while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line =~s/^\s+//;
    my @arr = split(/\s+/,$line);
    die "cannot recognize line: should be [rupture t, x, y, z, mag, T00, T01, T02, T10, T11, T12, T20, T21, T22, M00, M01, M02, M10, M11, M12, M20, M21, M22]\n" if($#arr != 23);
    my $fct = $arr[5] * $factor;
    print "$arr[6] $arr[7] $arr[8] $arr[9] $arr[10] $arr[11] $arr[12] $arr[13] $arr[14] $arr[2] $arr[3] $arr[4] $fct\n";
}
