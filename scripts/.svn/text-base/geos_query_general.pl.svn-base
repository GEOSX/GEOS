#!/usr/bin/perl
my $narg = $#ARGV + 1;

die "usage: $0 <session file> <query>\n" if($narg != 2);

#query = "Average Value";
#query = "Variable Sum";

#e.g., /Users/scottjohnson/Documents/workspace/gpac/test/full_tests/AdvectionSim
my $sname = shift;
my $query = shift;

print "import sys
RestoreSession(\"$sname\",0)
nplots = TimeSliderGetNStates();
start=0
end = nplots
for i in range(start, end):
    SetTimeSliderState(i)
    Query(\"$query\")
    print \"\%g\" \% GetQueryOutputValue()
sys.exit(\"\")
";
