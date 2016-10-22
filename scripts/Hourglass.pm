package Hourglass;
use strict;
use warnings;
use Math::MatrixReal;

require Exporter;
our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all'=> [qw(parse_himidlo get_his_los
				 slice_brute  hourglass_brute
				 compute_poly hourglass_poly
                                 quartic quadratic
				 random_images
                                 cov2ell
)]);
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw();



use Math::MatrixReal;
use Math::Complex;
use Math::Polynomial::Solve;

my $pi = atan(1,1)*4;

sub parse_himidlo {
  my $truthx = shift;
  my $truthy = shift;
  my $truthz = shift;

  my $hsh = {};
  for (@_) {
    chomp;
    my @ary = split /,/;
    my $iid = shift @ary;
    for my $key (qw(xhi yhi zhi  xmd ymd zmd  xlo ylo zlo)) {
      $hsh->{$iid}->{$key} = shift @ary;
      if ($key =~ /x/)    { $hsh->{$iid}->{$key} -= $truthx }
      if ($key =~ /y/)    { $hsh->{$iid}->{$key} -= $truthy }
      #if ($key =~ /[xy]/) { $hsh->{$iid}->{$key} *= 100000 }
      #if ($key eq 'zhi')  { $hsh->{$iid}->{$key}  =   10.0 }
      #if ($key eq 'zlo')  { $hsh->{$iid}->{$key}  =  -10.0 }
      if ($key =~ /z/)    { $hsh->{$iid}->{$key} -= $truthz }
    }
  }
  return $hsh;
}

sub get_his_los {
  my $hsh = shift;
  my @xhis=(); my @yhis=(); my @xlos=(); my @ylos=();
  my ($zhi, $zlo);
  for my $iid (@_) {
    my $h = $hsh->{$iid};
    push @xhis, $h->{xhi};
    push @yhis, $h->{yhi};
    push @xlos, $h->{xlo};
    push @ylos, $h->{ylo};
    $zhi = $h->{zhi};
    $zlo = $h->{zlo};
  }
  return (\@xhis, \@yhis, \@xlos, \@ylos, $zhi, $zlo);
}

sub slice_brute {
  my $xhis = shift; # these are array references
  my $yhis = shift;
  my $xlos = shift;
  my $ylos = shift;
  my $lam  = shift; # this is a scalar, 0=lo,1=hi

  my ($sx,$sy,$sxx,$syy,$sxy)=(0,0,0,0,0);
  my $n = @$xhis + 0; # all arrays should be the same length
  for my $i (0..$n-1) {
    my $x = $lam*$$xhis[$i] + (1-$lam)*$$xlos[$i];
    my $y = $lam*$$yhis[$i] + (1-$lam)*$$ylos[$i];
    $sx += $x;  $sxx += $x*$x;
    $sy += $y;  $syy += $y*$y;
    $sxy += $x*$y;
  }
  my $avgx = $sx/$n;
  my $avgy = $sy/$n;
  my $varx = ($sxx - $sx*$sx/$n) / ($n-1);
  my $vary = ($syy - $sy*$sy/$n) / ($n-1);
  my $cvar = ($sxy - $sx*$sy/$n) / ($n-1);
  my $m = Math::MatrixReal->new_from_rows( [ [$varx, $cvar],
                                             [$cvar, $vary] ]);
  my $e = $m->sym_eigenvalues();
  my $area = $pi * sqrt($e->element(1,1)) * sqrt($e->element(2,1));
  return ($avgx, $avgy, $area);
}

sub hourglass_brute {
  my $hsh = shift;
  my $nplanes = shift;
  my ($xhis, $yhis, $xlos, $ylos, $zhi, $zlo) = get_his_los($hsh, @_);

  my @axyzls;
  my @lams = ();
  for my $i (0..$nplanes) { push @lams, $i/$nplanes }
  for my $lam (@lams) {
    my ($avgx, $avgy, $area) = slice_brute($xhis, $yhis, $xlos, $ylos, $lam);
    my $zstr = sprintf "%.1f", $lam*$zhi + (1-$lam)*$zlo;
    push @axyzls, join(',', $area, $avgx, $avgy, $zstr, $lam);
  }

  my @saxyzls = sort {$a<=>$b} @axyzls;
  return $saxyzls[0];
}




sub compute_poly {
  my $xhis = shift; # these are array references
  my $yhis = shift;
  my $xlos = shift;
  my $ylos = shift;
  my $n = @$xhis+0;

  # compute all the sums and sums of cross products
  my %sum=(); # hash of sums
  for my $i (0..$n-1) {
    my %h = (); # hash these current values
    $h{xhi} = $$xhis[$i]; $sum{xhi} += $h{xhi};
    $h{yhi} = $$yhis[$i]; $sum{yhi} += $h{yhi};
    $h{xlo} = $$xlos[$i]; $sum{xlo} += $h{xlo};
    $h{ylo} = $$ylos[$i]; $sum{ylo} += $h{ylo};
    for my $ki (keys %h) { # tbd stop computing redundant lower diagonal
    for my $kj (keys %h) {
      $sum{$ki.$kj} += $h{$ki} * $h{$kj};
    }}
  }

  # now that we have all the sums, compute the covariances
  my %sig;
  for my $k1k2 (qw(xhixhi     xhixlo    xloxlo
                   yhiyhi     yhiylo    yloylo
                   xhiyhi xhiylo xloyhi xloylo)) {
    $k1k2 =~ /([xy][hl][io])([xy][hl][io])/;
    my $k1 = $1;
    my $k2 = $2;
    my $sumx  = $sum{$k1};
    my $sumy  = $sum{$k2};
    my $sumxy = $sum{$k1k2};
    $sig{$k1k2} = ($sumxy - $sumx*$sumy/$n) / ($n-1);
  }

  # var(x^i(\lambda)) = ax\lambda^s + bx\lambda + cx
  my $ax = $sig{xhixhi} - 2*$sig{xhixlo} +   $sig{xloxlo};
  my $bx =                2*$sig{xhixlo} - 2*$sig{xloxlo};
  my $cx =                                   $sig{xloxlo};

  # var(y^i(\lambda)) = ay\lambda^s + by\lambda + cy
  my $ay = $sig{yhiyhi} - 2*$sig{yhiylo} +   $sig{yloylo};
  my $by =                2*$sig{yhiylo} - 2*$sig{yloylo};
  my $cy =                                   $sig{yloylo};

  # covar(x,y) = axy\lambda^s + bxy\lambda + cxy
  my $axy= $sig{xhiyhi} - $sig{xhiylo} - $sig{xloyhi} +   $sig{xloylo};
  my $bxy=                $sig{xhiylo} + $sig{xloyhi} - 2*$sig{xloylo};
  my $cxy=                                                $sig{xloylo};

  # Now that we have quadratics(lambda) for varx/vary/cvar, we can
  # find the quartic coefficents of d(lambda) = varx*vary - cvar*cvar
  my $d4 = $ax*$ay - $axy*$axy;
  my $d3 = $ax*$by + $bx*$ay - 2*$axy*$bxy;
  my $d2 = $ax*$cy + $bx*$by + $cx*$ay - 2*$axy*$cxy - $bxy*$bxy;
  my $d1 = $bx*$cy + $cx*$by - 2*$bxy*$cxy;
  my $d0 = $cx*$cy - $cxy*$cxy;

  if ($ENV{DEBUG_HOURGLASS_POLY}) {
    # debug by comparing area for a specific lambda with the same
    # area computed by slice_brute
    my $l = 0.75;
    my $varx = $ax*$l*$l + $bx*$l + $cx;
    my $vary = $ay*$l*$l + $by*$l + $cy;
    my $cvar = $axy*$l*$l + $bxy*$l + $cxy;
    my $d_lam = $d4*$l*$l*$l*$l + $d3*$l*$l*$l + $d2*$l*$l + $d1*$l + $d0;
    my $polyarea = $pi * sqrt($d_lam);
    my ($avgx,$avgy,$area) = slice_brute($xhis, $yhis, $xlos, $ylos, $l);
    die "Polynomial calculation doesn't match brute force!!"
      unless (abs($polyarea-$area) < 1.0e-8);
  }

  return ($d4, $d3, $d2, $d1, $d0,
          $ax,  $bx,  $cx,
          $ay,  $by,  $cy,
          $axy, $bxy, $cxy);
}

sub quadratic {
  my $l = shift;
  return $_[0]*$l*$l + $_[1]*$l + $_[2];
}

sub quartic {
  my $l = shift;
  my $l2 = $l*$l;
  return $_[0]*$l2*$l2 + $_[1]*$l2*$l + $_[2]*$l2 + $_[3]*$l + $_[4];
}

sub hourglass_poly {
  my $hsh = shift @_;
  my $n = @_;
  my ($xhis, $yhis, $xlos, $ylos, $zhi, $zlo) = get_his_los($hsh, @_);

  # get the coefficients of the quartic polynomial d(lambda)
  my @dq = compute_poly($xhis, $yhis, $xlos, $ylos);
  my @d   = (@dq)[ 0..4];
  my @qx  = (@dq)[ 5..7];
  my @qy  = (@dq)[ 8..10];
  my @qxy = (@dq)[11..13];

  if ($d[0] == 0) {
    my $stophere=1;
  }

  # minimize by solving the cubic determinant
  my @roots = Math::Polynomial::Solve::cubic_roots(
	4*$d[0], 3*$d[1], 2*$d[2], $d[3]);
  my $lam = $roots[0]; # first root always real
  my $mal = 1-$lam;
  # other roots are complex or real together
  # make sure we don't have 3 real
  my $ambig = ''; # if 1 root, no ambiguity
  if (!ref $roots[1] && !ref $roots[2]) {
    my $stophere=1;
    # capture degenerate situation for later analysis
    local $" = ',';
    my @extrema;
    push @extrema, $pi*sqrt(quartic($roots[0], @d));
    push @extrema, $pi*sqrt(quartic($roots[1], @d));
    push @extrema, $pi*sqrt(quartic($roots[2], @d));
    if ($ENV{DEBUG_HOURGLASS_POLY}) {
      for my $iid (@_) { print "IID,$iid\n" }
      print "XHI,@$xhis\nYHI,@$yhis\nXLO,@$xlos\nYLO,@$ylos\n";
      print "ROOTS,@roots\n";
      print "VALUS,@extrema\n";
    }
    # find smallest extrema = absolute minimum
    $lam = $roots[0];
    my $ex  = $extrema[0];
    if ($extrema[1] < $ex) { $lam=$roots[1]; $ex=$extrema[1] }
    if ($extrema[2] < $ex) { $lam=$roots[2]; $ex=$extrema[2] }
    $mal = 1-$lam;
    my @sextrema = sort {$a<=>$b} @extrema;
    $ambig = $sextrema[1] - $sextrema[0];
    my $sambig = "$ambig";
    if ($sambig =~ /i/) {
      my $stophere=1;
    }
  }

  # OK now we have the optimal lambda. Interpolate the intersections sliced by
  # the plane at the optimal height
  my $sx=0; my $sy=0;
  for my $i (0..$n-1) {
    $sx += $lam*$xhis->[$i] + $mal*$xlos->[$i];
    $sy += $lam*$yhis->[$i] + $mal*$ylos->[$i];
  }
  my $avgx = $sx/$n;
  my $avgy = $sy/$n;
  my $z = $lam*$zhi + $mal*$zlo;
  my $dh = sqrt($avgx*$avgx, $avgy*$avgy);

  # Compute error-related terms
  my $varx = quadratic($lam, @qx);
  my $vary = quadratic($lam, @qy);
  my $cvar = quadratic($lam, @qxy);
  $varx /= $n*$n;
  $vary /= $n*$n;
  $cvar /= $n*$n;
  my $area = $pi*sqrt(quartic($lam, @d));
  # my $check = $pi*sqrt($varx*$vary - $cvar*$cvar); # should be == $area

  #my $inside_am = ($dh < $area/$n/$MAGIC_NUMBER);
  my $covar = Math::MatrixReal->new_from_rows( [ [$varx, $cvar],
						 [$cvar, $vary] ]);
  my ($icov, $ellipse_stat, $inside_e39, $inside_e90, $inside_e95)=(0)x5;
  my ($icov00, $icov11, $icov01) = (0)x3;
  if ($varx > 0 && $vary > 0) {
    $icov = $covar->inverse();
    $icov00 = $icov->element(1,1);
    $icov11 = $icov->element(2,2);
    $icov01 = $icov->element(1,2);
    $ellipse_stat = ($icov00 * $avgx*$avgx
                   + $icov11 * $avgy*$avgy
                 + 2*$icov01 * $avgx*$avgy);
    $inside_e39 = ($ellipse_stat < 1.000);
    $inside_e90 = ($ellipse_stat < 4.605);
    $inside_e95 = ($ellipse_stat < 5.991);

  # if ($n >= 50) {
  # $total_count++;
  # $total_inside += $inside_e90;
  # $total_avg = $total_inside / $total_count;
  # print STDERR "INSIDE_AM\t$n\t$total_inside\t$total_count\t$total_avg\n";
    # }
  }


  return (join ',', 'POLY', $n, $area, $avgx, $avgy, $dh, $z,
	            $lam, $varx, $vary, $cvar, $icov00, $icov11, $icov01,
	            $ellipse_stat, $inside_e39, $inside_e90, $inside_e95, $ambig, "\n");
}





sub random_images {
  my $n = shift;
  my %unused = ();
  my %used = ();
  for my $iid (@_) { $unused{$iid} = 1 }
  while (keys %used < $n) {
    my $nleft = keys %unused;
    my $i = rand()*$nleft;
    my $iid = (keys %unused)[$i];
    delete $unused{$iid}; # so we don't sample it again
    $used{$iid}++;
  }
  return keys %used;
}


# Convert 2D covariance to ellipse (1-sigma)
# Returns array of 7 elements in this order:
# major_radius, major_axis_x, major_axis_y (x/y are unit)
# minor_radius, minor_axis_x, minor_axis_y (x/y unit, perp. to maj)
# angle (of major axis, radians ccw from x axis)
sub cov2ell {
  my $varx = shift;
  my $vary = shift;
  my $cvar = shift;

  # compute eigenvecs/vals for 2x2
  my $matstr = "[ $varx $cvar ]\n[ $cvar $vary ]\n";
  my $cov = Math::MatrixReal->new_from_string($matstr);
  my ($evals, $evecs) = $cov->sym_diagonalize();
  my ($mj,$mn) = (1,2); # Major is first eval/vec
  # or it could be the other way around
  if ($evals->element(1,1) < $evals->element(2,1)) { $mj=2; $mn=1 }
  my $angle = atan2($evecs->element(2,$mj), $evecs->element(1,$mj));

  return (sqrt($evals->element($mj,1)), $evecs->element(1,$mj), $evecs->element(2,$mj),
	  sqrt($evals->element($mn,1)), $evecs->element(1,$mn), $evecs->element(2,$mn),
	  $angle);
}



1;
