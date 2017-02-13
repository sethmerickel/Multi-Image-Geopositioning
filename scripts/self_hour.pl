#! /usr/bin/perl
use FindBin;
use lib $FindBin::Bin;
use Hourglass ':all';

use Getopt::Std;
%opt = (n=>100, m=>25, k=>100, r=>100);
getopts('n:k:m:', \%opt);


select STDOUT; $| = 1; # autoflush
#$nplanes = 200; # plus 1!

# near 36N 117.5W 1700 in UTM11N
$truthx =  455000;
$truthy = 3984000;
$truthz =    1700;

@lines = (<>); # slurp up projector log
$hsh = parse_projector(0,0,0, @lines);

@all_iids = sort keys %$hsh;

print (join ',', qw(Niter Nbig Nsma Krpt HX HY HZ M_VX VY VZ CXY CXZ CYZ
                                                FPC_VX VY VZ CXY CXZ CYZ
                                                  N_VX VY VZ CXY CXZ CYZ
                                         MX MY MZ   VX VY VZ CXY CXZ CYZ), "\n");

for $rpt (1..$opt{r}) {
  @iids = random_images($opt{n}, @all_iids);
  $sumx  = $sumy  = $sumz  = 0;
  $sumxx = $sumyy = $sumzz = 0;
  $sumxy = $sumxy = $sumyz = 0;
  for $nnn (1..$opt{k}) {
    @some_iids = random_images($opt{m}, @iids);
    $hour = hourglass_poly($hsh, @some_iids);
    ($x,$y,$z) = (split /,/, $hour)[3,4,6];
    $x -= $truthx;
    $y -= $truthy;
    $z -= $truthz;
    $sumx += $x; $sumxx += $x*$x; $sumxy += $x*$y;
    $sumy += $y; $sumyy += $y*$y; $sumyz += $y*$z;
    $sumz += $z; $sumzz += $z*$z; $sumxz += $x*$z;
    #print "$nnn,$x,$y,$z\n";
  }
  $nnn = $opt{k};
  if ($nnn > 1) {
    $covxx = ($sumxx - $sumx*$sumx/$nnn) / ($nnn-1);
    $covyy = ($sumyy - $sumy*$sumy/$nnn) / ($nnn-1);
    $covzz = ($sumzz - $sumz*$sumz/$nnn) / ($nnn-1);
    $covxy = ($sumxy - $sumx*$sumy/$nnn) / ($nnn-1);
    $covxz = ($sumxz - $sumx*$sumz/$nnn) / ($nnn-1);
    $covyz = ($sumyz - $sumy*$sumz/$nnn) / ($nnn-1);
    $cov = mkmat(3,3, $covxx, $covxy, $covxz,
                      $covxy, $covyy, $covyz,
	              $covxz, $covyz, $covzz);
  }
  @hc = flatten($cov); # raw sample covariance of the k m-samples
  # apply finite population correction factor, no sqrt because in variance
  # space, not stdev space
  $cov *= ($opt{n}-1) / ($opt{n}-$opt{m});
  @fc = flatten($cov);
  # apply 1/sqrt(n) projection from m-samples to n-population, except not sqrt
  $cov *= $opt{m} / $opt{n};
  @nc = flatten($cov);

  $hour = hourglass_poly($hsh, @iids);
  ($x,$y,$z) = (split /,/, $hour)[3,4,6];
  $hourx = $x - $truthx;
  $houry = $y - $truthy;
  $hourz = $z - $truthz;

  $gp0 = mkmat(3,1, 1,1,1);
  ($gp, $cov, $refvar) = wvmig($hsh, $gp0, @iids);
  $migx = $gp->element(1,1) - $truthx;
  $migy = $gp->element(2,1) - $truthy;
  $migz = $gp->element(3,1) - $truthz;
  @mc = flatten($cov);

  print (join ',', $rpt, $opt{n}, $opt{m}, $opt{k},
	 $hourx, $houry, $hourz, $hc[0], $hc[4], $hc[8], $hc[1], $hc[2], $hc[5],
                                 $fc[0], $fc[4], $fc[8], $fc[1], $fc[2], $fc[5],
                                 $nc[0], $nc[4], $nc[8], $nc[1], $nc[2], $nc[5],
	 $migx,  $migy,  $migz,  $mc[0], $mc[4], $mc[8], $mc[1], $mc[2], $mc[5], "\n");

}


