#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0  <filename(s)>\n" if ( $narg < 1 );

my @filenames = ();

#PRINT HEADING
{
	my $filenamestr = "";
	for ( my $i = 0 ; $i < $narg ; $i++ ) {
		$filenames[$i] = shift;
		$filenamestr .= " " . $filenames[$i];
	}
	my $now = localtime time;
	print "*HEADING
Combined file generated from $filenamestr at $now
**
********************************** P A R T S **********************************
*PART, NAME=Part-Default
**
********************************** N O D E S **********************************
*NODE, NSET=ALLNODES
";

}

my @nodeoffsets = (0);

my $last_offset = 0;

#print nodes and record offsets
for ( my $i = 0 ; $i <= $#filenames ; $i++ ) {
	my $file = $filenames[$i];
	die "cannot open $file\n" if ( !open( IN, "<$file" ) );
	my $node_on = 0;
	my $nmax    = 0;
	{
		while (<IN>) {
			if (/^\*/) {
				if (/NODE/) {
					$node_on = 1;
				}
				else {
					$node_on = 0;
				}
			}
			else {
				if ( $node_on == 1 ) {
					my $line = $_;
					chomp($line);
					$line =~ s/^s+//;
					my @arr = split( /\,\s*/, $line );
					my $nodenumber = $arr[0] + $last_offset;
					print "$nodenumber, $arr[1], $arr[2], $arr[3]\n";
					if ( $nodenumber > $nmax ) {
						$nmax = $nodenumber;
					}
				}
			}
		}
	}
	close IN;
	$nodeoffsets[ $i + 1 ] = $nmax;
	$last_offset = $nmax;
}

print "**** WARNING - MATERIAL    1 UNDEFINED FOR THE FOLLOWING ELEMENTS ****
*ORIENTATION,NAME=SOR1,DEFINITION=COORDINATES
0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00
0,0.000000E+00
";

#determine where the element regions are located
my %elregions = ();
for ( my $i = 0 ; $i <= $#filenames ; $i++ ) {
	my $file = $filenames[$i];
	open( IN, "<$file" );
	my $el_on = 0;
	my $nmax  = 0;
	{
		while (<IN>) {
			if (/^\*/) {
				if (/^\*\s*ELEMENT/) {
					if (/ELSET\=(\S+)/) {
						my $elregion = $1;
						if ( exists( $elregions{$elregion} ) ) {
							push( @{ $elregions{$elregion} }, $i );
						}
						else {
							$elregions{$elregion} = [$i];
						}
						#print STDERR "$i contains $elregion\n";
					}
				}
			}
		}
	}
}

#print elements
$last_offset = 0;
my $elregion;
foreach $elregion ( sort( keys %elregions ) )
{
	my $fileindices = $elregions{$elregion};
	for ( my $i = 0 ; $i <= $#{$fileindices} ; $i++ ) {
		my $file = $filenames[ $$fileindices[$i] ];
		die "cannot open $file" if(!open( IN, "<$file" ));
		my $el_on = 0;
		my $nmax  = 0;
		{
			while (<IN>) {
				if (/^\*/) {
					if (/^\*\s*ELEMENT/) {
						my $line = $_;
						if (/ELSET\=(\S+)/) {
							if ( $1 eq $elregion ) {
								$el_on = 1;
								if ( $i == 0 ) {
									print $line;
								}
							}
							else {
								$el_on = 0;
							}
						}
					}
					else {
						$el_on = 0;
					}
				}
				else {
					if ( $el_on == 1 ) {
						my $line = $_;
						chomp($line);
						$line =~ s/^s+//;
						my @arr = split( /\,\s*/, $line );
						my $elnumber = $arr[0] + $last_offset;
						print "$elnumber";
						for ( my $j = 1 ; $j <= $#arr ; $j++ ) {
							my $nodenumber = $arr[$j] + $nodeoffsets[$$fileindices[$i]];
							print ", $nodenumber";
						}
						print "\n";
						if ( $elnumber > $nmax ) {
							$nmax = $elnumber;
						}
					}
				}
			}
		}
		close IN;
		$last_offset = $nmax;
	}
}

#print nodesets
my %nsets = ();
for ( my $i = 0 ; $i <= $#filenames ; $i++ ) {
	my $file = $filenames[$i];
	open( IN, "<$file" );
	my $on = 0;
	my $nsetcurr = "";
	{
		while (<IN>) {
			if (/^\*/) {
				if (/^\*\s*NSET\s*,\s*NSET\s*=\s*(\S+)/) {
					$nsetcurr = $1;
					$nsets{$nsetcurr} = [];
					$on = 1;
				}
				else {
					$on = 0;
				}
			}
			else {
				if ( $on == 1 ) {
					my $line = $_;
					chomp($line);
					$line =~ s/^s+//;
					my @arr = split( /\,\s*/, $line );
					for ( my $j = 0 ; $j <= $#arr ; $j++ ) {
						my $nodenumber = $arr[$j] + $nodeoffsets[$i];
						push(@{$nsets{$nsetcurr}},$nodenumber);
					}
				}
			}
		}
	}
	close IN;
}

my @keys = keys %nsets;
my @values = values %nsets;
while(@keys)
{
    my $key = pop(@keys);
    my $arr = pop(@values);
    print "*NSET,NSET=$key\n";
    for(my $i = 0; $i <= $#$arr; $i++)
    {
	print "$$arr[$i],\n";
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

