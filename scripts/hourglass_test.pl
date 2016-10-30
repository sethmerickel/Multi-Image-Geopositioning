#! /usr/bin/perl
use FindBin;
use lib $FindBin::Bin;
use Hourglass ':all';

select STDOUT; $| = 1; # autoflush


sub VERIFY {
  my $v0 = shift;
  my $v1 = shift;
  my $tol = shift;
  my $msg = "@_";
  if (abs($v0-$v1)<$tol) {
    printf "%-40s PASS\n", $msg;
  } else {
    printf "%-40s FAIL:\n", $msg;
    print "Difference is $v0 - $v1 = ".($v0-$v1)." > tol=$tol\n";
  }
}


$pi = atan2(1,1)*4;


for $case (0..2) {
  if ($case == 0) {
    @lines = qw(0,1,1,1,0,0,0,-1,-1,-1
		1,1,-1,1,0,0,0,-1,1,-1
		2,-1,-1,1,0,0,0,1,1,-1
		3,-1,1,1,0,0,0,1,-1,-1);
    $hsh = parse_himidlo(0,0,0, @lines);
    @iids = (0..3);
  } elsif ($case == 1) {
    # taken from *_ERROR0.sup from himidlo_utm11n.csv
    @lines = qw(0,454941.6,3984071.4,1710,454938.6,3984066.3,1700,454935.6,3984061.2,1690
		1,454939.2,3984063.6,1710,454937.6,3984065.3,1700,454936.1,3984067.0,1690
		2,454930.8,3984066.9,1710,454933.7,3984060.7,1700,454936.5,3984054.5,1690
		3,454929.2,3984063.4,1710,454934.8,3984068.6,1700,454940.5,3984073.9,1690
		4,454939.5,3984070.4,1710,454939.7,3984061.9,1700,454940.0,3984053.4,1690
		5,454932.9,3984065.4,1710,454935.2,3984066.1,1700,454937.5,3984066.7,1690
		6,454929.2,3984078.8,1710,454932.6,3984070.3,1700,454935.9,3984061.8,1690
		7,454932.4,3984063.3,1710,454938.8,3984067.1,1700,454945.1,3984070.9,1690
		8,454939.9,3984065.2,1710,454936.4,3984058.5,1700,454933.0,3984051.8,1690
		9,454932.7,3984063.7,1710,454930.8,3984064.1,1700,454928.8,3984064.5,1690);
    $hsh = parse_himidlo(454936.0104, 3984064.0306, 1700, @lines);
    @iids = (0..9);
  } elsif ($case == 2) { # degenerate bimodal case
      @lines = qw(0,0,0,1,0,0,0,2,2,-1
		  1,0,1,1,0,0,0,2,-1,-1
		  2,1,1,1,0,0,0,-1,-1,-1
                  3,1,0,1,0,0,0,-1,2,-1
                  4,-2,-2,1,0,0,0,0,0,-1
                  5,-2,1,1,0,0,0,0,-1,-1
                  6,1,1,1,0,0,0,-1,-1,-1
                  7,1,-2,1,0,0,0,-1,0,-1);
    $hsh = parse_himidlo(0,0,0, @lines);
    @iids = (0..7);
  } else {
    die "Case $case not defined!"
  }

  ($xhis, $yhis, $xlos, $ylos, $zhi, $zlo) = get_his_los($hsh, @iids);
  @dqs = compute_poly($xhis, $yhis, $xlos, $ylos);
  @d = (@dqs)[0..4];

  $l = 0.1;
  $d1 = $d[0]*$l*$l*$l*$l + $d[1]*$l*$l*$l + $d[2]*$l*$l + $d[3]*$l + $d[4];
  $d2 = quartic($l, @d);
  VERIFY($d1, $d2, 1e-12, "Case $case: Quartic computation");
  $a1 = $pi * sqrt($d1);

  $vx = quadratic($l, @dqs[5,6,7]);
  $vy = quadratic($l, @dqs[8,9,10]);
  $cv = quadratic($l, @dqs[11,12,13]);
  VERIFY($d1, $vx*$vy-$cv*$cv, 1.0e-12, "d vs vx*vy-cv*cv");

  ($avgx, $avgy, $a2) = slice_brute($xhis, $yhis, $xlos, $ylos, $l);

  VERIFY($a1, $a2, 1e-12, "Case $case: Area: Brute vs Poly");

  if ($case == 0) {
    # special tests for idealized case
    $d0 = quartic(0.5,@d);
    $a0 = $pi*sqrt($d0);
    VERIFY(0,$a0,1e-6,"0 ideal area (poly)");

    ($avgx, $avgy, $area) = slice_brute($xhis, $yhis, $xlos, $ylos, 0.5);
    VERIFY(0,$avgx,1e-12,"0 ideal xpos (brute)");
    VERIFY(0,$avgy,1e-12,"0 ideal ypos (brute)");
    VERIFY(0,$area,1e-12,"0 ideal area (brute)");

    $a1 = $pi * sqrt(quartic(1.0, @d));
    $vx = quadratic(1.0, @dqs[5,6,7]);
    $vy = quadratic(1.0, @dqs[8,9,10]);
    $cv = quadratic(1.0, @dqs[11,12,13]);
    VERIFY(4/3, $vx, 1.0e-12, "ideal xvar");
    VERIFY(4/3, $vy, 1.0e-12, "ideal yvar");
    VERIFY( 0,  $cv, 1.0e-12, "ideal cov");
    $a2 = $pi * sqrt($vx*$vy - $cv*$cv);
    VERIFY($a1, $a2, 1.0e-12, "ideal area: poly vs vars");
  }

  if ($case == 2) {
    # special test for degenerate bimodal

    # verify that 0123 and 4567 are both ideal
    $mode0 = hourglass_poly($hsh, (0,1,2,3));
    $mode1 = hourglass_poly($hsh, (4,5,6,7));
    ($a0,$z0) = (split ',', $mode0)[2,6];
    ($a1,$z1) = (split ',', $mode1)[2,6];
    VERIFY(0, $a0, 1.0e-6, "0 area bimode 1");
    VERIFY(0, $a1, 1.0e-6, "0 area bimode 2");
    VERIFY( 0.5, $z0, 1.0e-4,  "0.5 z bimode 1");
    VERIFY(-0.5, $z1, 1.0e-4, "-0.5 z bimode 2");


    # now verify that there are two identical minima
    ##### 2016-10-12 this is not working out -- test case
    ##### not properly balanced?
    $dmin0 = quartic(.25, @d);
    $dmin1 = quartic(.75, @d);
    VERIFY($dmin0, $dmin1, 1.0e-6, "symmetric bimodal");
    $bimode = hourglass_poly($hsh, (0..7));

    if ($ENV{PRINT_DEGENERATE}) {
      for ($l=0; $l<=1.005; $l+=0.01) {
	$dp = quartic($l, @d);
	$ap = $pi * sqrt($dp);
	$ab = slice_brute($xhis, $yhis, $xlos, $ylos, $l);
	VERIFY($ab, $ap, 1.0e-8, "Degen lam=$l");
	print (join ',', 'DEGENERATE',$l, $dp, "\n");
      }

    }

  }
}


@hmlns=qw(0,7157.081,30403.840,454938.607921,3984066.35492,1700,0.222058830365,-1.56634851079,0.73322290075,1.66965540612,-0.00876410952405,-0.495686854338,-1.74283,-0.0116751,-1.08811e-006,1.74052,-2.65088e-005,0.0224549,3.8686,0.025915,2.49187e-006,-3.86322,5.90554e-005,-0.0498419,-8.58577,-0.0504134,7.65224e-006,8.57508,-0.000131809,0.110626,-4.01949,-993811,-994964,-6662.4,-12832.7,-81.2924,9.14212,2.20584e+006,2.20856e+006,14799,28482.3,185.478,-27.6173,-4.89624e+006,-4.90172e+006,-29271.7,-63218,-412.223,
1,4003.374,32369.571,454937.669623,3984065.32994,1700,-0.0659355507611,-1.93276777754,-0.323798025788,1.93083854722,-0.0118683398978,-0.299146781349,-1.9596,0.00374246,1.17185e-006,1.95247,-3.57227e-005,0.0290517,4.8648,-0.00929085,-2.90745e-006,-4.84678,8.8676e-005,-0.0721174,-12.0771,0.023065,7.21715e-006,12.0316,-0.000220125,0.179023,4.42336,-993865,-997273,1900.09,-14834.4,36.7923,-10.9803,2.46716e+006,2.47578e+006,-4717.07,36824.8,-91.3319,27.2602,-6.12443e+006,-6.14626e+006,11710.4,-91413.4,226.723,
2,31565.161,29004.471,454933.718163,3984060.76787,1700,1.6495390705,0.209728866707,0.341869927292,0.0193504950112,1.43537910877,-0.88146764147,-1.69656,-0.0182686,-1.60107e-006,1.68366,2.63104e-005,0.0193784,-0.481183,-0.00518138,-4.54454e-007,0.477762,7.46626e-006,0.00549889,-0.136474,-0.00146956,-1.28909e-007,0.135571,2.11861e-006,0.00156038,-0.688722,-993774,-1.00126e+006,-10782,-11501.2,-120.575,-0.195599,-281997,-283954,-3057.73,-3263.61,-34.2148,-0.0554379,-80020.4,-80528.3,-867.162,-926.086,-9.70883,
3,23502.414,28728.281,454934.884394,3984068.67744,1700,1.44824223858,-0.357413226709,0.635704953734,0.000618050655248,1.43821279237,0.752613293019,-1.62073,0.020813,2.98127e-006,1.62229,1.16045e-005,0.0182229,0.629276,-0.00808117,-1.14517e-006,-0.629663,-4.63529e-006,-0.00707292,-0.244318,0.0031483,4.91776e-006,0.244374,1.80429e-006,0.0027451,-11.2152,-993772,-992693,12760.4,-11131.9,149.865,4.33003,385715,385460,-4954.95,4320.77,-58.3643,-4.41976,-149697,-149667,1930.48,-1676.96,22.6494,
4,3808.469,4477.120,454939.789771,3984061.92977,1700,-0.0487190980173,-1.20195819644,1.01963148757,1.56322340429,-0.010777656575,0.0425509236615,-1.57657,0.0125034,1.75365e-007,1.56265,-2.04831e-005,-0.0206127,2.19053,-0.0173725,-2.43058e-007,-2.17098,2.84567e-005,0.0286371,-3.04358,0.0241378,3.34378e-007,3.01614,-3.95336e-005,-0.0397854,1.42988,-993817,-1.00249e+006,7949.28,13232,-99.0735,-1.98671,1.3807e+006,1.39292e+006,-11045.2,-18383.2,137.642,2.75947,-1.9182e+006,-1.93541e+006,15346.9,25539.6,-191.225,
5,4057.641,1711.627,454935.26498,3984066.12525,1700,0.0198264297112,-1.95163825841,-0.12138307305,1.90535980992,-0.0159159344295,0.433702325086,-1.95426,-0.00822934,-8.30387e-007,1.95273,-3.55422e-005,-0.031193,2.67478,0.0112634,1.1363e-006,-2.67237,4.86399e-005,0.0426885,-3.66096,-0.0154162,-1.55502e-006,3.65722,-6.65674e-005,-0.0584205,1.60653,-993897,-994421,-4189.05,15890,72.6493,-2.19952,1.36018e+006,1.36106e+006,5733.53,-21745.9,-99.4234,3.01011,-1.86144e+006,-1.86288e+006,-7847.46,29759.9,136.064,
6,50134.019,4661.519,454932.618707,3984070.31502,1700,-0.315371612865,-1.18577513635,0.901714174845,1.44115986381,-0.00396311319564,0.48530954843,-1.52226,0.00743107,3.02541e-006,1.51972,-2.1832e-005,-0.0197643,-3.11206,0.0151919,6.18324e-006,3.10708,-4.46359e-005,-0.0404082,-6.36223,0.0310579,1.26493e-005,6.35242,-9.12555e-005,-0.0826146,-0.147509,-993812,-995299,4860.9,12947.1,-57.3603,-0.301672,-2.03185e+006,-2.03473e+006,9937.32,26470.3,-117.274,-0.617848,-4.15412e+006,-4.15966e+006,20315.2,54118.7,-239.766,
7,33272.483,5443.128,454938.815591,3984067.16805,1700,0.278376984006,-1.55339902444,-0.416240486877,1.38201055552,-0.00146315222519,0.881517799298,-1.63117,-0.0118636,-2.53654e-006,1.63846,-2.66044e-005,-0.020019,-1.46028,-0.0106208,-2.27196e-006,1.46703,-2.38204e-005,-0.0179245,-1.3073,-0.00950813,-2.03361e-006,1.31354,-2.13281e-005,-0.0160491,-8.46149,-993793,-989223,-7187.39,12093,96.2927,-7.5762,-889816,-885561,-6434.21,10827.8,86.2181,-6.78349,-796718,-792762,-5759.96,9694.96,77.1976,
8,21355.099,29706.420,454936.464544,3984058.54799,1700,0.291430166069,-1.39723208193,0.830236549294,1.53841351609,0.0157119350752,-0.541849804508,-1.65051,0.0147206,-2.00521e-006,1.63006,-2.29404e-005,0.0198861,1.3173,-0.0117488,1.6015e-006,-1.30121,1.83119e-005,-0.0158742,-1.05137,0.00937698,-1.27878e-006,1.0387,-1.46188e-005,0.0126717,4.07856,-993795,-1.00611e+006,8967.47,-12269,117.493,-3.25538,793305,802968,-7156.85,9793.83,-93.7896,2.59872,-633262,-640842,5711.82,-7817.99,74.8682,
9,37690.821,32753.270,454930.833913,3984064.15079,1700,-0.0253227729525,-1.98313891724,-0.0692151852301,1.93017816635,-0.0092348524453,-0.376989034053,-1.98324,0.00222976,9.77132e-007,1.96523,-3.64586e-005,0.0300004,2.92872,-0.00329276,-1.44186e-006,-2.90246,5.38442e-005,-0.0443077,-4.32495,0.00486254,2.11671e-006,4.28664,-7.95197e-005,0.0654381,2.41841,-993876,-1.00275e+006,1124.78,-15300.2,22.9548,-3.57076,1.46786e+006,1.48079e+006,-1661,22596.9,-33.9009,5.27962,-2.16788e+006,-2.18674e+006,2452.86,-33373.3,50.0708,);

$stop_at_hourmig=1; # for searching for a breakpoint
$ph = parse_partials(454936.0104, 3984064.0306, 1700, @hmlns);

# shift gp by ~1gsd, expect movement of ~1pix
$gp = $ph->{0}->{gp}->clone();
$gp->assign(1,1, $gp->element(1,1)+0.5);
$ip = wvg2i($ph, 0, $gp);
$ip -= $ph->{0}->{ip};
$dist = sqrt( (~$ip * $ip)->element(1,1) );
VERIFY(1, $dist, 0.2, "wvg2i GSD~0.5m");

# verify g2i with residual reported in hourmig.txt, NIMG=10, iid=0
$truth = mkmat(3,1, 454936.0104,   3984064.0306,   1700);
$gp    = mkmat(3,1, 454937.066757, 3984064.446951, 1699.945589);
$gp   -= $truth;
$proj  = wvg2i($ph, 0, $gp);
$res   = $ph->{0}->{ip} - $proj; # meas - proj
VERIFY(-2.60, $res->element(1,1), 0.01, "wvgwi res line");
VERIFY( 2.53, $res->element(2,1), 0.01, "wvgwi res samp");

# use the same gp0 as reported in hourmig NIMG=10
$gp0 = mkmat(3,1, 454937.717163, 3984065.544832, 1698.332910);
$gp0 -= $truth;
($gp1, $cov1, $refvar1) = wvmig($ph, $gp0, (0..9));
$cpp   = mkmat(3,1, 454937.066757, 3984064.446951,1699.945589);
$cpp -= $truth; # gp1 should match this
$err  = $gp1 - $cpp;
$dist = sqrt( (~$err*$err)->element(1,1));
VERIFY(0, $dist, 0.001, "wvmig matches cpp");

# send in gp that was output by cpp
($gp2, $cov2, $refvar2) = wvmig($ph, $cpp, (0..9));

($gp2, $cov2, $refvar2) = wvmig($ph, $gp1, (0..9));
$gp2 -= $gp1;
$dist = sqrt( (~$gp2 * $gp2)->element(1,1) );
VERIFY(0, $dist, 0.0001, "wvg2i 2nd iteration~0");





