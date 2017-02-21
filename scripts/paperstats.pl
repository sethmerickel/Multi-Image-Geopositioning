#! /bin/perl

use Hourglass('mkmat');

sub avg3 {
  my $sx=0;
  my $sy=0;
  my $sz=0;
  my $i=0;
  while ($i < @_) {
    $sx += $_[$i++];
    $sy += $_[$i++];
    $sz += $_[$i++];
  }
  my $k = @_/3;
  return ($sx/$k, $sy/$k, $sz/$k);
}

sub cele90 {
  my $c33 = shift;
  my $c22 = mkmat(2,2, $c33->element(1,1), $c33->element(1,2),
                       $c33->element(2,1), $c33->element(2,2));
  my $evals = $c22->sym_eigenvalues();
  my $e1 = $evals->element(1,1);
  my $e2 = $evals->element(2,1);
  my $emaj = sqrt($e1>$e2?$e1:$e2);
  my $emin = sqrt($e1<$e2?$e1:$e2);
  my $rho  = $emin / $emaj;
  my @CE90_coeffs = (+1.6439733846,
		     +0.0620678846,
		     -0.1849111597,
		     +1.1233898601,
		     -0.4982196970);
  my $poly_val = pop @CE90_coeffs; # quartic coefficient
  for (my $k=3; $k>=0; --$k) {
    $poly_val = $poly_val*$rho + (pop @CE90_coeffs);
  }
  my $ce90 = $poly_val * $emaj;
  my $le90 = sqrt($c33->element(3,3)) * 1.644874;
  return ($ce90, $le90);
}

$zero33 = mkmat(3,3, (0)x9);
$cov = $zero33->clone();

while (<>) {
  next if /ALG/;
  @ary = split /,/;
  $n = $ary[1];
  ($migx, $migy, $migz) = (@ary)[2,3,4];
  $refvar = $ary[5];
  ($vxx, $vxy, $vxz, $vyy, $vyz, $vzz) = (@ary)[6..11];
  ($plyx, $plyy, $plyz) = (@ary)[15,16,18];
  $migx -= 455000;  $plyx -= 455000;
  $migy -= 3984000; $plyy -= 3984000;
  $migz -= 1700;    $plyz -= 1700;
  print (join ',', 'ALL', $n, $refvar, $migx, $migy, $migz, $plyx, $plyy, $plyz, "\n");
  push @migs, $migx, $migy, $migz;
  #push @plys, $plyx, $plyy, $plyz;
  push @dhs, sqrt($migx*$migx + $migy*$migy);
  push @dzs,  abs($migz);
  $thiscov = mkmat(3,3, $vxx, $vxy, $vxz,
                        $vxy, $vyy, $vyz,
                        $vxz, $vyz, $vzz);
  $cov += $thiscov;

  if (@dhs==100 || $n == 1000) {
    $k = @dhs;
    @sdhs = sort {$a<=>$b} @dhs;
    @sdzs = sort {$a<=>$b} @dzs;
    $meas_ce90 = $sdhs[$k*0.9];
    $meas_le90 = $sdzs[$k*0.9];
    $cov *= 1.0/$k;
    ($pred_ce90,$pred_le90) = cele90($cov);
    ($adx, $ady, $adz) = avg3(@migs);
    print (join ',', 'AGG', $n, 
	   $adx, $ady, $adz,
	   $pred_ce90, $pred_le90,
	   $meas_ce90, $meas_le90, "\n");
    @dhs = @dzs = @migs = ();
    $cov = $zero33->clone();
  }
}
