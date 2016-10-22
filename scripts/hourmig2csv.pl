#! /usr/bin/perl
use FindBin;
use lib $FindBin::Bin;
use Hourglass ':all';

# 36N 117.5W 1700 in UTM11N
$truthx =  454936.0104;
$truthy = 3984064.0306;
$truthz = 1700;

print (join ',', qw(N DX DY DZ VX VY VZ REFV ELL), "\n");

while (<>) {
  if (/NIMG/)   { $n                    = (split)[-1]      }
  if (/UTM/)    { ($x,$y,$z)            = (split)[2,3,4]   }
  if (/COV/)    { ($varx, $vary, $cvar) = (split)[2, 5, 7] }
  if (/REFVAR/) { $refvar               = (split)[-1]      }

  if (/REFVAR/) {
    $dx = $x-$truthx;
    $dy = $y-$truthy;
    $dz = $z-$truthz;
    $dh = sqrt($dx*$dx + $dy*$dy);
    $d  = sqrt($dx*$dx + $dy*$dy + $dz*$dz);
    @ellipse = cov2ell($varx, $vary, $cvar);
    print (join ',', $n, $dx, $dy, $dz, $varx, $vary, $cvar, $refvar,
	   $dh/sqrt($n), $d/sqrt($n),
	   "\n");

  }
}
