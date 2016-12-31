#! /usr/bin/perl

use Math::MatrixReal;

use Getopt::Std;
%opt = ();
getopts('sf', \%opt);

#https://www.asprs.org/a/publications/proceedings/Sacramento2012/files/Dolloff.pdf
#interpolating not so finely here, just between (.75,1.9)-->(1.0,2.15)
sub compute_ce90 {
  my $m3x3 = shift;
  my $m2x2 = Math::MatrixReal->new(2,2);
  for my $r (1..2) {
  for my $c (1..2) {
    $m2x2->assign($r,$c, $m3x3->element($r,$c));
  }}
  my ($evals, $evecs) = $m2x2->sym_diagonalize();
  my $esma = $evals->element(1,1);
  my $ebig = $evals->element(2,1);
  my $r = $esma/$ebig;
  my $interp = ($r-0.75)*4; # scale 0.75-1 to 0-1
  my $k = $interp*2.15 + (1-$interp)*1.9;
  return $k * sqrt($ebig);
}


for $f (@ARGV) { # output files from migpoly.pl
  open MIGPOLY, $f;
  while (<MIGPOLY>) {
    next if /ALG,N/; # skip the header line
    @ary = split /,/;
    shift @ary; # 'MIG'
    $n = shift @ary; # N
    for $key (qw(DX DY DZ REF ECXX ECXY ECXZ ECYY ECYZ ECZZ)) {
      $val = shift @ary;
      if ($key eq 'DX') { $val -= 455000;  $dx = $val }
      if ($key eq 'DY') { $val -= 3984000; $dy = $val }
      if ($key eq 'DZ') { $val -= 1700;    $dz = $val }
      $stat->{$f}->{$n}->{$key} += $val;
    }
    $stat->{$f}->{$n}->{D3} += sqrt($dx*$dx+$dy*$dy+$dz*$dz);
    $stat->{$f}->{$n}->{DH} += sqrt($dx*$dx+$dy*$dy        );
    $stat->{$f}->{$n}->{DV} += sqrt(                $dz*$dz);
    $stat->{$f}->{$n}->{DXDX} += $dx*$dx;
    $stat->{$f}->{$n}->{DXDY} += $dx*$dy;
    $stat->{$f}->{$n}->{DXDZ} += $dx*$dz;
    $stat->{$f}->{$n}->{DYDY} += $dy*$dy;
    $stat->{$f}->{$n}->{DYDZ} += $dy*$dz;
    $stat->{$f}->{$n}->{DZDZ} += $dz*$dz;
    $stat->{$f}->{$n}->{num} += 1;
  }
  close MIGPOLY;
}

@hary1 = ('');
for $f (@ARGV) { push @hary1, $f; push @hary1, ('')x30; }
print (join ',', @hary1, "\n");

@hary2 = ('N');
for $f (@ARGV) { push @hary2, qw(DX DY DZ D3 DH DV REF
                                 ECXX ECXY ECXZ ECYY ECYZ ECZZ
                                 ACXX ACXY ACXZ ACYY ACYZ ACZZ
		       CE90 LE90 ERH1 ERH2 ERZ EVX EVY EVZ
				 ARH1 ARH2 ARZ AVX AVY AVZ) }
print (join ',', @hary2, "\n");

@ns = sort {$a<=>$b} keys %{$stat->{$ARGV[0]}};
for $n (@ns) {
  @pary = ($n);
  for $f (@ARGV) {
    $s = $stat->{$f}->{$n};
    $num = $s->{num};
    next if $num<=1;
    $fact = 1;
    if ($opt{s}) { $fact *= sqrt($n); }
    if ($opt{f}) { $fact *= 1.0/sqrt((1000-$n)/(1000-1)); }
    for $key (qw(DX DY DZ D3 DH DV REF)) {
      push @pary, $s->{$key}/$num * $fact;
    }
    for $key (qw(ECXX ECXY ECXZ ECYY ECYZ ECZZ)) {
      push @pary, $s->{$key}/$num;
    }

    ($exx, $exy, $exz, $eyy, $eyz, $ezz) = (@pary)[-6,-5,-4,-3,-2,-1];

    $axx = ($s->{DXDX} - $s->{DX}*$s->{DX}/$num) / ($num-1);
    $axy = ($s->{DXDY} - $s->{DX}*$s->{DY}/$num) / ($num-1);
    $axz = ($s->{DXDZ} - $s->{DX}*$s->{DZ}/$num) / ($num-1);
    $ayy = ($s->{DYDY} - $s->{DY}*$s->{DY}/$num) / ($num-1);
    $ayz = ($s->{DYDZ} - $s->{DY}*$s->{DZ}/$num) / ($num-1);
    $azz = ($s->{DZDZ} - $s->{DZ}*$s->{DZ}/$num) / ($num-1);
    push @pary, $axx, $axy, $axz, $ayy, $ayz, $azz;

    $cov = Math::MatrixReal->new(3,3);
    $cov->assign(1,1,$exx); $cov->assign(1,2,$exy); $cov->assign(1,3,$exz);
    $cov->assign(2,1,$exy); $cov->assign(2,2,$eyy); $cov->assign(2,3,$eyz);
    $cov->assign(3,1,$exz); $cov->assign(3,2,$eyz); $cov->assign(3,3,$ezz);
    $ce90 = compute_ce90($cov);
    $le90 = sqrt($cov->element(3,3))*1.64;
    push @pary, $ce90*$fact, $le90*$fact;

    ($evals, $evecs) = $cov->sym_diagonalize();
    push @pary, sqrt($evals->element(1,1)) * $fact;
    push @pary, sqrt($evals->element(2,1)) * $fact;
    push @pary, sqrt($evals->element(3,1)) * $fact;
    push @pary, $evecs->element(1,3);
    push @pary, $evecs->element(2,3);
    push @pary, $evecs->element(3,3);

    $cov->assign(1,1,$axx); $cov->assign(1,2,$axy); $cov->assign(1,3,$axz);
    $cov->assign(2,1,$axy); $cov->assign(2,2,$ayy); $cov->assign(2,3,$ayz);
    $cov->assign(3,1,$axz); $cov->assign(3,2,$ayz); $cov->assign(3,3,$azz);
    ($evals, $evecs) = $cov->sym_diagonalize();
    push @pary, sqrt($evals->element(1,1)) * $fact;
    push @pary, sqrt($evals->element(2,1)) * $fact;
    push @pary, sqrt($evals->element(3,1)) * $fact;
    push @pary, $evecs->element(1,3);
    push @pary, $evecs->element(2,3);
    push @pary, $evecs->element(3,3);
  }
  next if (@pary == 1);
  print (join ',', @pary, "\n");

  #last if $n >= 200;
}
