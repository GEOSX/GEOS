#!/usr/bin/perl -w

#==========================================================
#
#  lune.pl
#  Perl script to write a shell script for GMT plotting of moment tensor source types.
#  Carl Tape, carltape@gi.alaska.edu, May 2012
#
#  Reference:
#  W. Tape and C. Tape, "A geometric setting for moment tensors," Geophysical J. International, 2012
#
#  Last tested 8-22-2012 with GMT 4.5.8
#  note: GMT's psmeca does not currently handle full moment tensor beachballs properly;
#        a custom version is needed to use this option
#  
#==========================================================

#use Math::Trig;
my $narg = $#ARGV + 1;
die "usage: $0 <label(s) space separated>\nNote: a file of the format \"dfiles/beachpts_%s_pts.dat\" where %s is the label should exist in the \"dfiles\" sub-directory of the current working directory\n" if($narg == 0);

my @ftags = ();
my @inds = ();
for(my $i = 0; $i < $narg; $i++)
{
    my $label = shift;
    push(@ftags, $label);
    push(@inds, $i);
    #"Ford2009","Foulger2004","Minson2007","Minson2008","Walter2009","Walter2010","Pesicek2012","Baig2010","Sileny2006","Sileny2008","Sileny2009","Dreger2012");
}

$cshfile = "lune.csh";
$bindir = "/opt/local/lib/gmt4/bin";
$fontno = "1"; 

open(CSH,">$cshfile");
print CSH "$bindir/gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.3c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 18 FRAME_PEN 2p TICK_PEN 2p\n";

$R = "-R-30/30/-90/90";
#$R = "-R-30/30/-90/90r";
$origin = "-X2 -Y1";
$xtick1 = 10; $ytick1 = 10;
$xtick2 = 5; $ytick2 = 5;
#$B = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":WesN";
$B = "-Ba${xtick1}f${xtick2}g${xtick1}:\" \":/a${ytick1}f${xtick2}g${ytick1}:\" \":wesn";

$wid = "3i";
#$J = "-JA0/0/$wid"; $title = "Lambert equal-area ($J)"; $ftag = "lambert";
#$J = "-JB0/0/0/90/$wid"; $title = "Albers equal-area ($J)"; $ftag = "albers";

#$J = "-JY0/$wid";  $title1 = "Cylindrical equal-area ($J)"; $ftag = "cylindrical";
#$J = "-JI0/$wid";  $title1 = "Sinusoidal equal-area ($J)"; $ftag = "sinusoidal";
#$J = "-JKf0/$wid"; $title1 = "Eckert IV equal-area ($J)"; $ftag = "eckert4";
#$J = "-JK0/$wid";  $title1 = "Eckert VI equal-area ($J)"; $ftag = "eckert6";
#$J = "-JW0/$wid";  $title1 = "Mollweide equal-area ($J)"; $ftag = "mollewide";
$J = "-JH0/$wid";   $title1 = "Hammer equal-area ($J)"; $ftag = "hammer";

# colors
#$magenta = "148/0/211";
$magenta = "160/32/240";
$orange = "255/165/0";
$red = "255/0/0";
$blue = "30/144/255";
$cyan = "0/255/255";
$green = "50/205/50";
$sienna = "160/82/45";
$brown = "139/69/16";
$black = "0/0/0";
$white = "255/255/255";

$lgray = 200;
$dgray = 120;

# KEY COMMAND
$iplot = 1;  # =0 (reference lune), =1 (dots from published studies), =2 (reference beachballs)
$lplot = 1;  # =1-2: reference MTs on the lune (iplot=2 only)
$kplot = 3;  # =1-4: orientation of MT at center of lune (iplot=2 only)
if($iplot==2) {
  $psfile = "lune_${ftag}_iplot${iplot}_lplot${lplot}_kplot${kplot}.ps";
} else {
  $psfile = "lune_${ftag}_iplot${iplot}.ps";
}
$ipatch = 1;  # three shaded patches on the lune
$icrack = 1;  # nu=0.25 arc between crack points
$ilam2 = 1;   # lam2 = 0 arc between dipoles
$ilegend = 0; # legend for data points
if ($iplot==0) {$plot_ref_points = 1; $plot_ref_labels = 1;}
if ($iplot==1) {$plot_ref_points = 0; $plot_ref_labels = 1;}
if ($iplot==2) {$plot_ref_points = 0; $plot_ref_labels = 0;}
if ($iplot==1) {$ilegend = 1;}

$clune = $lgray;
if($ipatch==0) {$clune = $sienna;}

print CSH "$bindir/psbasemap $J $R $B -G$clune -K -V -P $origin > $psfile\n"; # START

$pdir = "./dfiles/";

# plot patches
if ($ipatch==1) {
  $fname = "$pdir/beach_patch_01.lonlat";
  print CSH "$bindir/psxy $fname -G$dgray -J -R -K -O -V >>$psfile\n";
  $fname = "$pdir/beach_patch_02.lonlat";
  print CSH "$bindir/psxy $fname -G255 -J -R -K -O -V >>$psfile\n";
  print CSH "$bindir/psbasemap $J $R $B -K -O -V >> $psfile\n";
}

# plot arcs
# dev, iso+DC, iso, iso, CDC nu=0.25, CDC nu=0
$lwid = 3;
$fname1 = "$pdir/beach_arc_01.lonlat";  # deviatoric (equator)
$fname2 = "$pdir/beach_arc_02.lonlat";  # iso+DC (center longitude)
$fname3 = "$pdir/beach_arc_03.lonlat";  # bottom of +ISO patch
$fname4 = "$pdir/beach_arc_04.lonlat";  # top of -ISO patch
$fname5 = "$pdir/beach_arc_05.lonlat";  # CDC nu=0.25
$fname6 = "$pdir/beach_arc_06.lonlat";  # CDC nu=0 (between dipoles)
if ($iplot==1) {
  $W = "-W2p,0,--";
  print CSH "$bindir/psxy $fname1 $W -J -R -K -O -V >>$psfile\n";
  print CSH "$bindir/psxy $fname2 $W -J -R -K -O -V >>$psfile\n";
  print CSH "$bindir/psxy $fname3 $W -J -R -K -O -V >>$psfile\n";
  print CSH "$bindir/psxy $fname4 $W -J -R -K -O -V >>$psfile\n";
  if($icrack==1) {print CSH "$bindir/psxy $fname5 $W -J -R -K -O -V >>$psfile\n";}
  if($ilam2==1) {print CSH "$bindir/psxy $fname6 $W -J -R -K -O -V >>$psfile\n";}
} else {
  if($ipatch==1) {@cols = ($magenta,$red,$blue,$blue,$black,$blue);}
  else           {@cols = ($magenta,$orange,$red,$white,$black,$blue);}
  if($ipatch==1) {print CSH "$bindir/psxy $fname1 -W${lwid}p,$cols[0] -J -R -K -O -V >>$psfile\n";}
  if($ipatch==1) {print CSH "$bindir/psxy $fname2 -W${lwid}p,$cols[1] -J -R -K -O -V >>$psfile\n";}
  print CSH "$bindir/psxy $fname3 -W${lwid}p,$cols[2] -J -R -K -O -V >>$psfile\n";
  print CSH "$bindir/psxy $fname4 -W${lwid}p,$cols[3] -J -R -K -O -V >>$psfile\n";
  if($icrack==1) {print CSH "$bindir/psxy $fname5 -W${lwid}p,$cols[4] -J -R -K -O -V >>$psfile\n";}
  if($ilam2==1) {print CSH "$bindir/psxy $fname6 -W${lwid}p,$cols[5] -J -R -K -O -V >>$psfile\n";}
}

# plot lune reference points and labels
if ($plot_ref_points || $plot_ref_labels) {
  $fname = "$pdir/beach_points.lonlat";
  if ($plot_ref_points) {
    $csize = 12;
    print CSH "$bindir/psxy $fname -N -Sc${csize}p -W1p,0/0/0 -G255 -J -R -K -O -V >>$psfile\n";
  }
  if ($plot_ref_labels) {
    $fsize = 14;
    $fontno = 1;
    open(IN,$fname); @plines = <IN>; close(IN);
    for ($i = 1; $i <= @plines; $i++) {
      ($plon,$plat,$plab,$Dx,$Dy) = split(" ",$plines[$i-1]);
      #print "\n--$plon -- $plat-- $plab --";
      $D = "-D${Dx}p/${Dy}p";
      print CSH "$bindir/pstext -N -J -R -K -O -V $D >>$psfile<<EOF\n$plon $plat $fsize 0 $fontno CM $plab\nEOF\n";
    }
    #print CSH "awk '{print \$1,\$2,$fsize,0,0,\"CM\",\$3}' $fname | pstext -N -J -R -K -O -V >> $psfile\n";
  }
}

if ($iplot==1) {
  # moment tensors from various studies
  $csize = 8;  # size of dots
  @csizes = ($csize,$csize,$csize,$csize,$csize,$csize,$csize,$csize/2,$csize,$csize,$csize,$csize/2);
  @cols = ($black,$orange,$red,$white,$green,$magenta,$black,$green,$red,$orange,$green,$cyan);
  #@ftags = ($label);#"Ford2009","Foulger2004","Minson2007","Minson2008","Walter2009","Walter2010","Pesicek2012","Baig2010","Sileny2006","Sileny2008","Sileny2009","Dreger2012");
  #@inds = (8,2,3,4,7);     # volcanic/geothermal + Baig
  #@inds = (8,9,10,11);    # induced
  #@inds = (0);#1,2,3,5,6);    # TapeTape2012 figure
  #@inds = (9..11);       
  #@inds = 12;

    for ($i = 1; $i <= @inds; $i++) {
      $j = $inds[$i-1];
      $cz = $csizes[$j-1];
      $fname = sprintf("$pdir/beachpts_%s_points.dat",$ftags[$j-1]);
      print CSH "$bindir/psxy $fname -N -Sc${cz}p -W0.5p,0/0/0 -G$cols[$j-1] -J -R -K -O -V >>$psfile\n";
    }

#   # Minson vs GCMT
#   @ftits = ("Minson2007 (n=14)","GCMT (n=14)");
#   @cols = ($red,$cyan);
#   @inds = (1,2);
#   $fname1 = "/home/carltape/papers/SOURCE_INVERSION/DATA/MinsonGCMT_Minson.dat";
#   $fname2 = "/home/carltape/papers/SOURCE_INVERSION/DATA/MinsonGCMT_GCMT.dat";
#   print CSH "awk '{print \$8,\$9}' $fname1 | psxy -N -Sc${csize}p -W0.5p,0/0/0 -G$cols[0] -J -R -K -O -V >> $psfile\n";
#   print CSH "awk '{print \$8,\$9}' $fname2 | psxy -N -Sc${csize}p -W0.5p,0/0/0 -G$cols[1] -J -R -K -O -V >> $psfile\n";

#    # Dreger et al. 2012
#    $cptfile = "color.cpt";
#    #print CSH "makecpt -Crainbow -T50/100/5 -D > $cptfile\n";
#    print CSH "makecpt -Cseis -T50/100/5 -D -I > $cptfile\n";
#    $ilegend = 0;
#    $fname = "/home/carltape/papers/SOURCE_INVERSION/DATA/bsldreger_fmt_lune.dat";
#    $Fmin = 40;  # try 40,70,90
#    `awk '\$3 > $Fmin' $fname > dtemp`;
#    $nplot = `wc dtemp | awk '{print \$1}'`; chomp($nplot);
#    print "\n$nplot FMTs with F > $Fmin\n";
#    print CSH "awk '{print \$1,\$2,\$3}' dtemp | psxy -N -Sc${csize}p -W0.5p,0/0/0 -C$cptfile -J -R -K -O -V >> $psfile\n";
#    $Dscale = "-D0/1/2/0.2";
#    $Bscale = "-B10f5:\"F-test significance\": -Eb10p";
#    print CSH "psscale -C$cptfile $Dscale $Bscale -Xa3.5 -Ya6 -V -K -O >> $psfile\n";

} elsif ($iplot==2) {
  # reference beachballs on the lune
  $cmtinfo = "-Sm0.5 -L0.5p/0/0/0 -G255/0/0 -N";
  $cmtfile = sprintf("$pdir/beachballs_ilune%i_iref%i_psmeca",$lplot,$kplot);
  print CSH "$bindir/psmeca $cmtfile $J $R $cmtinfo -K -O -V >> $psfile\n";
} 

#-----------------------------

$J_title = "-JX1i";  # -JM7i
$R_title = "-R0/1/0/1";
$otitle1 = "-Xa3.5 -Ya6";

# legend for plotting published studies
if($iplot==1 && $ilegend==1) {
  $x0 = 0; $y0 = 1.2; $dy = 0.3;
  for ($i = 1; $i <= @inds; $i++) {
    $j = $inds[$i-1];
    $cz = $csizes[$j-1];
    $x = $x0 + 0.2;
    $y = $y0 - ($i-1)*$dy;
    # number of points (for legend)
    $fname = sprintf("$pdir/beachpts_%s_points.dat",$ftags[$j-1]);
    $nplot = `wc $fname | awk '{print \$1}'`; chomp($nplot);
    $lab = sprintf("%s (n=%i)",$ftags[$j-1],$nplot);
    print CSH "$bindir/psxy -N -Sc${cz}p -W0.5p,0/0/0 -G$cols[$j-1] $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n$x0 $y\nEOF\n";
    print CSH "$bindir/pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n $x $y 12 0 $fontno LM $lab\nEOF\n";
  }

#$x = -30; $y = rad2deg(asin(1/sqrt(3)));
#print CSH "psxy -N -Sc${csize}p -W1p,0/0/0 -G255/165/0 -J -R -K -O -V >>$psfile<<EOF\n$x $y\nEOF\n";
}

#-----------------------------

# optional: plot a title
$otitle1 = "-Xa-1 -Ya9.0";
$otitle2 = "-Xa-1 -Ya8.7";
if (0==1) {
  $title1 = "Representation of source types on the fundamental lune";
  $title2 = "(W. Tape and C. Tape, 2012, GJI, \"A geometric setting for moment tensors\")";
  #$title1 = "Berkeley Seismological Laboratory full moment tensor catalog (n = $nplot, Fsig > $Fmin)";
  #$title2 = "Dreger, Chiang, Ford, Walter, 2012, Monitoring Research Review";
  print CSH "$bindir/pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0 0 14 0 $fontno LM $title1\nEOF\n";
  print CSH "$bindir/pstext -N $R_title $J_title $otitle2 -O -V >>$psfile<<EOF\n 0 0 11 0 $fontno LM $title2\nEOF\n";
} else {
  # any command with no -K should work (this one will not plot anything)
  print CSH "$bindir/pstext $R_title $J_title $otitle2 -O -V >>$psfile<<EOF\n -1 -1 11 0 $fontno LM TEST\nEOF\n"; 
}

close (CSH);
system("csh -f $cshfile");

#system("$bindir/ps2pdf $psfile");
# you may need to install gv to view (or use something else)
#system("gv $psfile &");

# to make composite pdf file:
# for file in `ls lune_hammer_iplot2_*.ps` ; do ps2pdf $file ; done ; pdcat -r lune_hammer_iplot2*pdf all_lune_hammer_iplot2.pdf
# pdcat -r lune_hammer_iplot0_kplot1.pdf lune_hammer_iplot1_kplot1.pdf all_lune_hammer_iplot2.pdf all_lune.pdf

#==================================================

