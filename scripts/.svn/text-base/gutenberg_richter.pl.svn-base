#!/opt/local/bin/perl
use strict;

my $narg = $#ARGV + 1;
if ( $narg == 0 ) {
	print
"usage: $0 <filename of time(years) magnitude tuples> <optional: time_period> <optional: number of bins to auto determine> <optional: bin start points>\n";
	exit();
}

my $filename = shift;
my $time_period = -1;
if($narg > 1)
{
    $time_period = shift;
}

#--------------------
# GET THE BINS
#--------------------
my @bins   = ();
my @counts = ();
if ( $narg > 3 ) 
{
    $narg-=2;
    
    #get the bins
    for(my $i = 0;$i < $narg; $i++)
    {
		my $val = shift;
		push(@bins,$val);
		push(@counts,0);
    }
}
else 
{
	my $binnum = 20;
	if ( $narg == 3 ) 
	{
		$binnum = shift;
	}
	$binnum++;
	my $min = 1e100;
	my $max = -1e100;
	die "cannot open file $filename" if ( !open( IN, "<$filename" ) );
	while (<IN>) 
	{
		my $line = $_;
		$line =~ s/^\s+//;
		chomp($line);
		my @arr = split(/\s+/,$line);
		die "incorrect number of fields in record: $line" if($#arr != 1);
		$line = $arr[1];
		if ( $line > $max ) 
		{
			$max = $line;
		}
		if ( $line < $min )
		{
			$min = $line;
			die "a0=$arr[0] a1=$arr[1] --> $_" if($line eq 'inf' or $line eq '-inf');
		}
	}
	close IN;
	my $binsize = ( $max - $min ) / ( $binnum - 1 );
	my $current_bin = $min;
	for ( my $i = 0 ; $i < $binnum ; $i++ ) 
	{
		push( @bins,   $current_bin );
		push( @counts, 0 );
		$current_bin += $binsize;
	}
	$bins[$#bins] = $max;
}

#--------------------
# READ THE FILE
#--------------------
my $count = 0;
my $tmax = -1;
die "cannot open file $filename" if ( !open( IN, "<$filename" ) );
while (<IN>) 
{
	my $line = $_;
	$line =~ s/^\s+//;
	chomp($line);
	my @arr = split(/\s+/,$line);
	die "incorrect number of fields in record: $line" if($#arr != 1);
	$line = $arr[1];
	if($arr[0] > $tmax)
	{
		$tmax = $arr[0];
	}
	for ( my $i = 0 ; $i <= $#bins ; $i++ )
	{
		if ( $line >= $bins[$i] ) 
		{
			$counts[$i]++;
		}
		else 
		{
			last;
		}
	}
}
close IN;

die "max time is invalid: $tmax" if($tmax <= 0);

if($time_period < 0)
{
    $time_period = $tmax;
}

#-----------------------------------
# CALCULATE STATISTICS AND PRINT OUTPUT
#-----------------------------------

my $nbins     = $#bins + 1;
my $nmin      = 4;
my $nmintotal = 2 * $nmin + 1;

if ( $nbins < $nmintotal || 1 ) {

	#single linear regression line
	my $xy_sum = 0;
	my $x_sum  = 0;
	my $x2_sum = 0;
	my $y_sum  = 0;
	for ( my $i = 0 ; $i <= $#bins ; $i++ ) {

		#calculate N per annum
		my $nnorm = $counts[$i] / ($time_period);
		die "cannot have counts <= 0 : $i -> $counts[$i] ($nbins and $bins[$i]) :: nnorm = $nnorm\n" if($nnorm <= 0);
		my $val   = log10($nnorm);
		print "$i $bins[$i] $counts[$i] $nnorm $val\n";
		$xy_sum += $val * $bins[$i];
		$x_sum  += $bins[$i];
		$y_sum  += $val;
		$x2_sum += $bins[$i] * $bins[$i];
	}
    my $inv_n = 1 / ($#bins + 1);
    my $b_value = -($xy_sum - ($inv_n * $x_sum * $y_sum)) / ($x2_sum - ($inv_n * $x_sum * $x_sum));
    my $a_value = $inv_n * ($y_sum + $b_value * $x_sum);
    print "b-value: $b_value a-value: $a_value\n"; 
}
else 
{
	#optimal bilinear regression
	my @rms12   = ( 0, 0 );
	my @nbins12 = ( 0, 0 );
	my @beta1s  = ( 0, 0 );

	my @xx_sum = ( 0, 0 );
	my @xy_sum = ( 0, 0 );
	my @x_mean = ( 0, 0 );
	my @y_mean = ( 0, 0 );

	my @x = ();
	my @y = ();

	#---------------------------
	#get first estimate
	$nbins12[0] = $nmin;
	$nbins12[1] = $nbins - $nbins12[0];
	for ( my $iset = 0 ; $iset < 2 ; $iset++ ) 
	{

		#set beginning and last+1 indices
		my $i0 = 0;
		my $i1 = $nbins12[0];
		if ( $iset > 0 ) {
			$i0 = $nbins12[0];
			$i1 = $nbins;
		}

		#get initial values of state variables
		for ( my $i = $i0 ; $i < $i1 ; $i++ ) {
			my $nnorm = $counts[$i] / ($time_period);
			my $val   = log10($nnorm);
			$x[$i] = $bins[$i];
			$y[$i] = $val;
			print "$bins[$i] $counts[$i] $nnorm $val\n";

			$xx_sum[$iset] += $x[$i] * $x[$i];
			$xy_sum[$iset] += $x[$i] * $y[$i];
			$x_mean[$iset] += $x[$i];
			$y_mean[$iset] += $y[$i];
		}
		$x_mean[$iset] /= $nbins12[$iset];
		$y_mean[$iset] /= $nbins12[$iset];

		#calculate initial rms
		my $ssxy =
		  &SSxy( $xy_sum[$iset], $x_mean[$iset], $y_mean[$iset],
			$nbins12[$iset] );
		my $ssxx = &SSxx( $xx_sum[$iset], $x_mean[$iset], $nbins12[$iset] );
		my $b1 = &Beta1( $ssxy, $ssxx );
		my $b0 = &Beta0( $b1, $x_mean[$iset], $y_mean[$iset] );
		$beta1s[$iset] = -$b1;

		#die "b-value: $beta1s[$iset] [$i0 , $i1 ] -> (xmean $x_mean[$iset] ymean $y_mean[$iset] xxsum $xx_sum[$iset] xysum $xy_sum[$iset]\n" if($iset == 1);
		for ( my $i = $i0 ; $i < $i1 ; $i++ ) {
			my $val = &Ytilde( $x[$i], $b0, $b1 );
			$val -= $y[$i];
			$rms12[$iset] += $val * $val;
		}
		$rms12[$iset] /= $nbins12[$iset] - 2;
		$rms12[$iset] = sqrt( $rms12[$iset] );
	}

	#---------------------------
	#get successive estimates
	for ( my $i = 0 ; $i < ( $nbins - $nmintotal ) ; $i++ ) {
		my $n1 = $i + $nmin + 1;
		my $n2 = $nbins - $n1;
		my $xa = $x[$n1];
		my $ya = $y[$n1];

		$xx_sum[0] += $xa * $xa;
		$xx_sum[1] -= $xa * $xa;

		$xy_sum[0] += $xa * $ya;
		$xy_sum[1] -= $xa * $ya;

		$x_mean[0] = ( ( $n1 - 1 ) * $x_mean[0] + $xa ) / $n1;
		$x_mean[1] = ( ( $n2 + 1 ) * $x_mean[1] - $xa ) / $n2;

		$y_mean[0] = ( ( $n1 - 1 ) * $y_mean[0] + $ya ) / $n1;
		$y_mean[1] = ( ( $n2 + 1 ) * $y_mean[1] - $ya ) / $n2;

		#calculate trms
		my @trms = ( 0, 0 );
		my @tnbins = ( $n1, $n2 );
		my @tb1 = ( 0, 0 );
		for ( my $iset = 0 ; $iset < 2 ; $iset++ ) 
		{
			#set beginning and last+1 indices
			my $i0 = 0;
			my $i1 = $tnbins[0];
			if ( $iset > 0 ) 
			{
				$i0 = $tnbins[0];
				$i1 = $nbins;
			}
			my $ssxy = &SSxy(
				$xy_sum[$iset], $x_mean[$iset],
				$y_mean[$iset], $tnbins[$iset]
			);
			my $ssxx = &SSxx( $xx_sum[$iset], $x_mean[$iset], $tnbins[$iset] );
			my $b1 = &Beta1( $ssxy, $ssxx );
			my $b0 = &Beta0( $b1, $x_mean[$iset], $y_mean[$iset] );
			$tb1[$iset] = -$b1;
			for ( my $i = $i0 ; $i < $i1 ; $i++ ) 
			{
				my $val = &Ytilde( $x[$i], $b0, $b1 );
				$val -= $y[$i];
				$trms[$iset] += $val * $val;
			}
			$trms[$iset] /= $tnbins[$iset] - 2;
			$trms[$iset] = sqrt( $trms[$iset] );
		}
		my $srms_curr      = $rms12[0] + $rms12[1];
		my $srms_candidate = $trms[0] + $trms[1];

  		#print "trying ... $i ($tnbins[0] $tnbins[1]) $srms_curr $srms_candidate -> ";
		if ( $srms_candidate < $srms_curr ) 
		{
			for ( my $iset = 0 ; $iset < 2 ; $iset++ ) 
			{
				$rms12[$iset]   = $trms[$iset];
				$beta1s[$iset]  = $tb1[$iset];
				$nbins12[$iset] = $tnbins[$iset];
			}

			#print "smallest!\n";
		}
		else 
		{

			#print "reject\n";
		}
	}
	print "b-values: $beta1s[0] $beta1s[1] $nbins12[0] $nbins12[1] $rms12[0] $rms12[1]\n";
}

sub SSxx {
	my $xx_sum = shift;
	my $x_mean = shift;
	my $n      = shift;
	my $ret    = $xx_sum - ( $n * $x_mean * $x_mean );
	return $ret;
}

sub SSxy {
	my $xy_sum = shift;
	my $x_mean = shift;
	my $y_mean = shift;
	my $n      = shift;
	my $ret    = $xy_sum - ( $n * $x_mean * $y_mean );
	return $ret;
}

sub Beta1 {
	my $SSxy = shift;
	my $SSxx = shift;
	my $ret  = $SSxy / $SSxx;
	return $ret;
}

sub Beta0 {
	my $beta1  = shift;
	my $x_mean = shift;
	my $y_mean = shift;
	my $ret    = $y_mean - ( $beta1 * $x_mean );
	return $ret;
}

sub Ytilde {
	my $x     = shift;
	my $beta0 = shift;
	my $beta1 = shift;
	my $ret   = $beta0 + ( $beta1 * $x );
	return $ret;
}

sub log10 {
	my $n = shift;
	return log($n) / log(10);
}
