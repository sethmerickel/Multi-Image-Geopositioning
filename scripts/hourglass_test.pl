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


for $case (0..1) {
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
}







