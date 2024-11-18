#!/usr/bin/perl
use strict;

my @arr = ();
{
    my $line = <STDIN>;
    $line =~s/^\s+//;
    chomp($line);
    @arr = split(/\s+/,$line);
}
die "usage: $0 <stdin of x y z Txx Txy ... Tzz>\n" if($#arr != 11);

my @T = ([$arr[3],$arr[4],$arr[5]],[$arr[6],$arr[7],$arr[8]],[$arr[9],$arr[10],$arr[11]]);
my @q = (rand(),rand(),rand(), rand());
{
    my $sum = 0;
    for(my $i = 0; $i <= $#q; $i++)
    {
	$sum += $q[$i]*$q[$i];
    }
    $sum = 1.0/sqrt($sum);
    for(my $i = 0; $i <= $#q; $i++)
    {
	$q[$i] *= $sum;
    }
}
my @A = &QuaternionToTensor(@q);
my @Ttilde = &RotateTensor(\@A, \@T);

print "$arr[0] $arr[1] $arr[2] ${$Ttilde[0]}[0] ${$Ttilde[0]}[1] ${$Ttilde[0]}[2] ${$Ttilde[1]}[0] ${$Ttilde[1]}[1] ${$Ttilde[1]}[2] ${$Ttilde[2]}[0] ${$Ttilde[2]}[1] ${$Ttilde[2]}[2]";

sub RotateTensor
{
    my $pA = shift;
    my $pT = shift;

    my @A = @$pA;
    my @T = @$pT;

    my @AT = ([0,0,0],[0,0,0],[0,0,0]);
    for(my $i = 0; $i < 3; $i++)
    {
	for(my $j = 0; $j < 3; $j++)
	{
	    for(my $k = 0; $k < 3; $k++)
	    {
		${$AT[$i]}[$j] += ${$A[$i]}[$k] * ${$T[$k]}[$j];
	    }
	}
    }
    return @AT;

    die "${$AT[0]}[0] ${$AT[0]}[1] ${$AT[0]}[2]\n";
    my @ATAt = ([0,0,0],[0,0,0],[0,0,0]);
    for(my $i = 0; $i < 3; $i++)
    {
	for(my $j = 0; $j < 3; $j++)
	{
	    for(my $k = 0; $k < 3; $k++)
	    {
		${$ATAt[$i]}[$j] += ${$AT[$i]}[$k] * ${$A[$j]}[$k];
	    }
	}
    }



    return @ATAt;    
}

sub QuaternionToTensor
{
    my $xx = $_[0]*$_[0] + $_[1]*$_[1] - $_[2]*$_[2] - $_[3]*$_[3];
    my $yy = $_[0]*$_[0] - $_[1]*$_[1] + $_[2]*$_[2] - $_[3]*$_[3];
    my $zz = $_[0]*$_[0] - $_[1]*$_[1] - $_[2]*$_[2] + $_[3]*$_[3];

    my $xy = 2 * ($_[1]*$_[2] - $_[0]*$_[3]);
    my $xz = 2 * ($_[1]*$_[3] + $_[0]*$_[2]);

    my $yx = 2 * ($_[1]*$_[2] + $_[0]*$_[3]);
    my $yz = 2 * ($_[2]*$_[3] - $_[0]*$_[1]);

    my $zx = 2 * ($_[1]*$_[3] - $_[0]*$_[2]);
    my $zy = 2 * ($_[2]*$_[3] + $_[0]*$_[1]);

    return ([$xx,$xy,$xz],[$yx,$yy,$yz],[$zx,$zy,$zz]);
}
