#!/usr/perl
use strict;

my $narg = $#ARGV + 1;
die "usage: $0 <filename>" if($narg != 1);

my $filename = shift;
die "cannot open $filename" if(!open(IN,"<$filename"));

my $p_on = 0;
die "nothing to read" if(eof(IN));

#FIRST, GET PSTR
my $pstr = "";
{
    my $line = <IN>;
    chomp($line);
    my $ok = 1;
    while($ok)
    {
	my @arr = split(/\<Parameter /,$line);
	#print "$arr[0]\n";
	if($#arr > 0)
	{
	    my @arr2 = split(/\/\>/,$arr[1]);
	    $pstr .= " $arr2[0]";
	    while($#arr2 == 0)
	    {
		die "reached end of file before end of parameter node!" if(eof(IN));
		$line = <IN>;
		chomp($line);
		@arr2 = split(/\/\>/,$line);
		$pstr .= " $arr2[0]";
	    }
	    $line = $arr[1];
	    chomp($line);
	}
	else
	{
	    if(eof(IN))
	    {
		$ok = 0;
	    }
	    else
	    {
		$line = <IN>;
		chomp($line);
	    }
	}
    }
}
close IN;

#SECOND, PRINT EVERYTHING ELSE
open(IN,"<$filename");
{
    my $found = 0;
    my $line = <IN>;
    chomp($line);
    my $ok = 1;
    while($ok)
    {
	my @arr = split(/\<Parameter /,$line);
	print "$arr[0]\n";
	if($#arr > 0)
	{
	    @arr = split(/\/\>/,$line);
	    while($#arr == 0)
	    {
		die "reached end of file before end of parameter node!" if(eof(IN));
		$line = <IN>;
		chomp($line);
		@arr = split(/\/\>/,$line);
	    }
	    if($found == 0)
	    {
		&ExtractParameters($pstr);
		$found = 1;
	    }
	    $line = $arr[1];
	    chomp($line);
	}
	else
	{
	    if(eof(IN))
	    {
		$ok = 0;
	    }
	    else
	    {
		$line = <IN>;
		chomp($line);
	    }
	}
    }
}
close IN;
exit;

sub ExtractParameters
{
    my $pstr = shift;
    my @arr = split(/\"\s+/,$pstr);
    print "\t<Parameters>\n";
    for(my $i = 0; $i <= $#arr; $i++)
    {
	my @arr2 = split(/\=/,$arr[$i]);
	my $name = $arr2[0];
	my $value = $arr2[1];
	$value =~s/\"//g;
	print "\t\t<Parameter name=\"$name\" value=\"$value\"/>\n";
    }
    print "\t</Parameters>\n";
}
