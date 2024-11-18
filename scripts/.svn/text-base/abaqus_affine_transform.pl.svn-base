#!/bin/perl
use strict;

my $narg = $#ARGV + 1;
die
"usage: $0 <file name> <translate-x = 0> <opt: translate-y = 0> <opt: translate-z = 0> <opt: Txx = 1> <opt: Txy = 0> <opt: Txz = 0> <opt: Tyx = 0> <opt: Tyy = 1> <opt: Tyz = 0> <opt: Tzx = 0> <opt: Tzy = 0> <opt: Tzz = 1>\n"
  if ( $narg < 1 || $narg > 13 );

#--------------------------------
#READ IN ARGUMENTS
#--------------------------------

my @dx = ( 0, 0, 0 );
my @T = ( 1, 0, 0, 0, 1, 0, 0, 0, 1 );

my $fname = shift;
die "cannot open $fname\n" if ( !open( IN, "<$fname" ) );

my $ii = 1;
for ( my $i = 0 ; $i < 3 ; $i++ ) {
	if ( $ii >= $narg ) {
		last;
	}
	else {
		$dx[$i] = shift;
		$ii++;
	}
}
for ( my $i = 0 ; $i < 9 ; $i++ ) {
	if ( $ii >= $narg ) {
		last;
	}
	else {
		$T[$i] = shift;
		$ii++;
	}
}

#--------------------------------
#READ ABAQUS
#--------------------------------
{
	my $node_on = 0;
	my $icurr   = 0;
	while (<IN>) {
		if (/^\s*\*\s*(\S+)/) {
			if (/NODE/) {
				$node_on = 1;
			}
			else {
				$node_on = 0;
			}
			print $_;
		}
		else {
			if ($node_on) {
				my $line = $_;
				chomp($line);
				$line =~ s/^\s+//;
				my @arr = split( /\,\s*/, $line );
				if ( $#arr != 3 ) {
					die "node_on but this doesn't look like a node: $_";
				}
				my @x = ( $arr[1], $arr[2], $arr[3] );
				$icurr++;

				my @xx = &Transform( \@x, \@T, \@dx );
				print "$arr[0], $xx[0], $xx[1], $xx[2]\n";
			}
			else {
				print $_;
			}
		}
	}
}
close IN;

#--------------------------------
#TRANSFORM
#--------------------------------

sub Transform {
	my $xptr  = shift;
	my $Tptr  = shift;
	my $dxptr = shift;

	my @x  = @{$xptr};
	my @T  = @{$Tptr};
	my @dx = @{$dxptr};

	my @xnew = ( $dx[0], $dx[1], $dx[2] );
	my $ii = 0;
	for ( my $i = 0 ; $i < 3 ; $i++ ) {
		for ( my $j = 0 ; $j < 3 ; $j++ ) {
			$xnew[$i] += $T[$ii] * $x[$j];
			$ii++;
		}
	}
	return @xnew;
}
