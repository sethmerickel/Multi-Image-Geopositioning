#! /usr/bin/perl

use FindBin;
use lib $FindBin::Bin;
use Hourglass ':all';

#https://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html
@thresh  = (1.0000, 6.2513, 7.8147, 11.3448);
@percent = (0.1988, 0.90,   0.95,    0.99);

$truth = mkmat(3,1, 455000, 3984000, 1700);

while (<>) {
  next if /ALG/;
  @ary = split /,/;
  #  0  1 2 3 4  5  6    7      8      9    10     11
  # ALG,N,X,Y,Z,REF,VARX,CVARXY,CVARXZ,VARY,CVARYZ,VARZ
  $nimg = $ary[1];
  last if $nimg == 1000;
  $xyz = mkmat(3,1, (@ary)[2,3,4]);
  $dxyz = $xyz - $truth;

  $dxyz *= 1.0/sqrt( (1000-$nimg)/(1000-1) );

  $cov = mkmat(3,3, (@ary)[6,7,8, 7,9,10, 8,10,11]);
  #$refvar = $ary[5];
  #if ($refvar > 1) {
  #  $cov *= 1.0/$refvar;
  #}
  $voc = $cov->inverse();
  $ell = ~$dxyz * $voc * $dxyz;
  $ellipse_constant = $ell->element(1,1);

  $n++;
  for $i (0..$#thresh) {
    if ($ellipse_constant < $thresh[$i]) {
      $count[$i] += 1;
      $ct->{$nimg}->{$i} += 1;
    }
  }
}

print "$n observations\n";
for $i (0..$#thresh) {
  $should = $percent[$i];
  $is = $count[$i] / $n;
  printf "should be %.2f less than %5.2f;   is: %.2f\n", $should, $thresh[$i], $is;
}

exit;

for $nimg (sort {$a<=>$b} keys %$ct) {
  @pary = ($nimg);
  for $i (0..$#thresh) {
    $should = $percent[$i];
    $is = $ct->{$nimg}->{$i} / 100;
    push @pary, $is;
  }
  print (join ',', @pary, "\n");
}
