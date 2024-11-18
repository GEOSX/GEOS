#!/usr/bin/perl
my $narg = $#ARGV + 1;
die "usage: $0 <filename>
   filename line format: T00, T01, T02, T10, T11, T12, T20, T21, T22, X, Y, Z
   where [X,Y,Z] is the source location and [Tij] are the components of the Eigenvectors of the moment tensor
" if($narg != 1);

my $filename = shift;
die "cannot open $filename\n" if(!open(IN,"<$filename"));

print "#include \"colors.inc\"

#declare S11 = intersection { sphere { <0,0,0>, 0.5 pigment{ color Black } } box { <0,0,-1>, <1,1,1> pigment { color Black} } }
#declare S00 = intersection { sphere { <0,0,0>, 0.5 pigment{ color Black } } box { <0,0,-1>, <-1,-1,1> pigment { color Black} } }
#declare S01 = intersection { sphere { <0,0,0>, 0.5 pigment{ color White } } box { <0,0,-1>, <1,-1,1> pigment { color White} } }
#declare S10 = intersection { sphere { <0,0,0>, 0.5 pigment{ color White } } box { <0,0,-1>, <-1,1,1> pigment { color White} } }

#declare BEACHBALL = union { object { S11 } object { S00 } object { S10 } object { S01 } }
";

while(<IN>)
{
    my $line = $_;
    chomp($line);
    $line =~s/^\s+//;
    my @arr = split(/\s+/,$line);
    die "cannot recognize line: should be [T, x, m]\n" if($#arr != 12);
    print " object { BEACHBALL scale <$arr[12], $arr[12], $arr[12]> matrix <$arr[0], $arr[1], $arr[2], $arr[3], $arr[4], $arr[5], $arr[6], $arr[7], $arr[8], $arr[9], $arr[10], $arr[11]> }
";
}
close IN;

