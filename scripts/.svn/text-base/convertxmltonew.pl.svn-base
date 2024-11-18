#/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <filename>" if ( $narg != 1 );

my $filename = shift;
die "cannot open $filename" if ( !open( IN, "<$filename" ) );

my %elementRegionAttributes = ();
$elementRegionAttributes{"name"} = 1;
$elementRegionAttributes{"elementtype"} = 1;
$elementRegionAttributes{"hgDamp"} = 1;
$elementRegionAttributes{"hgStiff"} = 1;
$elementRegionAttributes{"abaqusID"} = 1;

my %contactAttributes = ();
$contactAttributes{"penetrationTol"} = 1;
$contactAttributes{"cosMinTol"} = 1;
$contactAttributes{"active"} = 1;

my $errd = 0;
while(<IN>)
{
	my $line = $_;
	$line=~s/setname=/setnames=/;
	if(/<ElementRegion /)
	{
		#die "ER found\n";
		my @arr = split(/<ElementRegion /);
		my $spc = $arr[0];
		my $heading = "ElementRegion";
		my %kvs = ();
		&ReadAttributes($heading, $line, \%kvs);
		print $spc;
		print "<ElementRegion";
		while(my($k,$v) = each %kvs)
		{
			#print "key: $k --> value: $v\n";
			my $ok = 0;
			if(exists $elementRegionAttributes{$k})
			{
				$ok = 1;
			}
			if($ok == 1)
			{
				print " $k=$v";
			}
		}
		print ">\n$spc$spc<LinearElasticMaterial";
		while(my($k,$v) = each %kvs)
		{
			my $ok = 1;
			if(exists $elementRegionAttributes{$k})
			{
				$ok = 0;
			}
			if($ok == 1)
			{
				print " $k=$v";
			}
		}
		print "/>\n$spc</ElementRegion>\n";
		$errd = 1;
	}
	elsif(/<Contact /)
	{
		my @arr = split(/<Contact /);
		my $spc = $arr[0];
		my $heading = "Contact";
		my %kvs = ();
		&ReadAttributes($heading, $line, \%kvs);
		print $spc;
		print "<Contact";
		while(my($k,$v) = each %kvs)
		{
			my $ok = 0;
			if(exists $contactAttributes{$k})
			{
				$ok = 1;
			}
			if($ok == 1)
			{
				print " $k=$v";
			}
		}
		print ">\n$spc$spc<PenaltyCoulombJoint";
		while(my($k,$v) = each %kvs)
		{
			my $ok = 1;
			if(exists $contactAttributes{$k})
			{
				$ok = 0;
			}
			if($ok == 1)
			{
				print " $k=$v";
			}
		}
		print "/>\n$spc</Contact>\n";		
	}	
	else
	{
		print $line;
	}
}

close IN;

sub ReadAttributes {
	my $heading = shift;
	my $line = shift;
	my $hashRefPtr = shift;

	my $p_on = 0;
	die "nothing to read" if ( eof(IN) );

	#FIRST, GET ATTRIBUTES ALL IN A SINGLE STRING
	my $pstr = "";
	{
		chomp($line);
		my $ok = 1;
		while ($ok) {
			my @arr = split( /\<$heading /, $line );
			if ( /\<$heading / ) {
				#you found the desired header
				$line = $arr[1];

				chomp($line);
				$line=~s/^\s+//;

				if($line=~m/\/>/)
				{
					my @arr2 = split( /\/>/, $line );
					my $str = $arr2[0];
					chomp($str);
					$str=~s/^\s+//;			
					$pstr .= " $str";
					$ok = 0;
				}
				else
				{
					$pstr.=" $line";
				}
				while ($ok) {
					die "reached end of file before end of parameter node!"
					  if ( eof(IN) );
					$line = <IN>;

					chomp($line);
					$line=~s/^\s+//;

					if($line=~m/\/>/)
					{
						my @arr2 = split( /\/>/, $line );
						my $str = $arr2[0];
						chomp($str);
						$str=~s/^\s+//;			
						$pstr .= " $str";
						$ok = 0;
					}
					else
					{
						$pstr.=" $line";
					}
				}
				last;
			}
			else {
				#keep looking for the header
				if ( eof(IN) ) {
					$ok = 0;
				}
				else {
					$line = <IN>;
				}
			}
		}
	}
	
	#SECOND, SPLIT ATTRIBUTES INTO KEY-VALUE PAIRS
	#die "$heading ::::::: $pstr\n";
	my @kvs = split("=",$pstr);
	if($#kvs >= 0)
	{
		my $keyValue = $kvs[0];
		$keyValue=~s/^\s+//;
		chomp($keyValue);		
		for(my $i = 1; $i <= $#kvs; $i++)
		{
			my $valueStr = $kvs[$i];
			my @arrSplit = split(/ /,$valueStr);
			if($#arrSplit >= 0)
			{
				my $valueValue = $arrSplit[0];
				for(my $j = 1; $j < $#arrSplit; $j++)
				{
					$valueValue .= " $arrSplit[$j]";
				}
				#print "key: $keyValue --> value: $valueValue\n";
				$valueValue=~s/^\s+//;
				chomp($valueValue);
				$$hashRefPtr{$keyValue} = $valueValue;
				$keyValue = $arrSplit[$#arrSplit];
				$keyValue=~s/^\s+//;
				chomp($keyValue);
			}
		}
	}
}
