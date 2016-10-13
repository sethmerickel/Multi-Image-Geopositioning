#! /usr/bin/perl
use FindBin;
use lib $FindBin::Bin;
use Hourglass ':all';

select STDOUT; $| = 1; # autoflush

if (@ARGV) {     # if you provide an input file with columns like the example below
  @lines = (<>); # slurp it in
  for (@lines) { push @iids, (split /,/)[0] }
} else {         # else use this as an example
  # taken from *_ERROR0.sup from himidlo_utm11n.csv
  # note xmid,ymid,zmid are ignored
  #         iid  xhi      yhi      zhi   xmid     ymid     zmid  xlo       ylo     zlo
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
  @iids = (0..9);
}

$hsh = parse_himidlo(454936.0104, 3984064.0306, 1700, @lines);
($xhis, $yhis, $xlos, $ylos, $zhi, $zlo) = get_his_los($hsh, @iids);
@dqs = compute_poly($xhis, $yhis, $xlos, $ylos);
@d   = (@dqs)[0..4];
@qx  = (@dqs)[5..7];
@qy  = (@dqs)[8..10];
@qxy = (@dqs)[11..13];

print (join ',', qw(LAM VARX VARY CVAR Z X0 Y0 X1 Y1), "\n");

for ($lam=0; $lam <= 1; $lam += 0.25) { # or smaller increments
  $mal = 1-$lam;
  $varx = quadratic($lam, @qx);
  $vary = quadratic($lam, @qy);
  $cvar = quadratic($lam, @qxy);
  $z = $lam*$zhi + $mal*$zlo;
  @xys = ();
  for $i (0..$#iids) {
    $x = $lam*$$xhis[$i] + $mal*$$xlos[$i];
    $y = $lam*$$yhis[$i] + $mal*$$ylos[$i];
    push @xys, $x, $y;
  }
  print (join ',', $lam, $varx, $vary, $cvar, $z, @xys, "\n");
}
