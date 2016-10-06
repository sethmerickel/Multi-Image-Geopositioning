#! /usr/bin/perl
use Math::MatrixReal;
use Math::Complex;
use Math::Polynomial::Solve;

sub slice_brute {
  my $xhis = shift; # these are array references
  my $yhis = shift;
  my $xlos = shift;
  my $ylos = shift;
  my $lam  = shift; # this is a scalar, 0=lo,1=hi

  my ($sx,$sy,$sxx,$syy,$sxy)=(0,0,0,0,0);
  my $n = @$xhis + 0; # all arrays should be the same length
  for $i (0..$n-1) {
    $x = $lam*$$xhis[$i] + (1-$lam)*$$xlos[$i];
    $y = $lam*$$yhis[$i] + (1-$lam)*$$ylos[$i];
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
  my @xhis=(), @yhis=(), @xlos=(), @ylos=();
  my $zhi, $zlo;
  for $iid (@_) {
    my $h = $hsh->{$iid};
    push @xhis, $h->{xhi};
    push @yhis, $h->{yhi};
    push @xlos, $h->{xlo};
    push @ylos, $h->{ylo};
    $zhi = $h->{zhi};
    $zlo = $h->{zlo};
  }
  my @axyzls;
  for my $lam (@lams) {
    my ($avgx, $avgy, $area) = slice_brute(\@xhis, \@yhis, \@xlos, \@ylos, $lam);
    $zstr = sprintf "%.1f", $lam*$zhi + (1-$lam)*$zlo;
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
  for $i (0..$n-1) {
    my %h = (); # hash these current values
    $h{xhi} = $$xhis[$i]; $sum{xhi} += $h{xhi};
    $h{yhi} = $$yhis[$i]; $sum{yhi} += $h{yhi};
    $h{xlo} = $$xlos[$i]; $sum{xlo} += $h{xlo};
    $h{ylo} = $$ylos[$i]; $sum{ylo} += $h{ylo};
    for $ki (keys %h) { # tbd stop computing redundant lower diagonal
    for $kj (keys %h) {
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
  my $n = @_;
  my @xhis=(), @yhis=(), @xlos=(), @ylos=();
  my $zhi, $zlo;
  for $iid (@_) {
    my $h = $hsh->{$iid};
    push @xhis, $h->{xhi};
    push @yhis, $h->{yhi};
    push @xlos, $h->{xlo};
    push @ylos, $h->{ylo};
    $zhi = $h->{zhi};
    $zlo = $h->{zlo};
  }

  # get the coefficients of the quartic polynomial d(lambda)
  my @dq = compute_poly(\@xhis, \@yhis, \@xlos, \@ylos);
  @d   = (@dq)[ 0..4];
  @qx  = (@dq)[ 5..7];
  @qy  = (@dq)[ 8..10];
  @qxy = (@dq)[11..13];

  # minimize by solving the cubic determinant
  my @roots = Math::Polynomial::Solve::cubic_roots(
	4*$d[0], 3*$d[1], 2*$d[2], $d[3]);
  my $lam = $roots[0]; # first root always real
  my $mal = 1-$lam;
  # other roots are complex or real together
  # make sure we don't have 3 real
  my $ambig = ''; # if 1 root, no ambiguity
  if (!ref $roots[1] && !ref $roots[2]) {
    $stophere=1;
    # capture degenerate situation for later analysis
    local $" = ',';
    my @extrema;
    push @extrema, $pi*sqrt(quartic($roots[0], @d));
    push @extrema, $pi*sqrt(quartic($roots[1], @d));
    push @extrema, $pi*sqrt(quartic($roots[2], @d));
    if ($ENV{DEBUG_HOURGLASS_POLY}) {
      for $iid (@_) { print "IID,$iid\n" }
      print "XHI,@xhis\nYHI,@yhis\nXLO,@xlos\nYLO,@ylos\n";
      print "ROOTS,@roots\n";
      print "VALUS,@extrema\n";
    }
    # find smallest extrema = absolute minimum
    $lam = $roots[0];
    $ex  = $extrema[0];
    if ($extrema[1] < $ex) { $lam=$roots[1]; $ex=$extrema[1] }
    if ($extrema[2] < $ex) { $lam=$roots[2]; $ex=$extrema[2] }
    $mal = 1-$lam;
    my @sextrema = sort {$a<=>$b} @extrema;
    $ambig = $sextrema[1] - $sextrema[0];
    $sambig = "$ambig";
    if ($sambig =~ /i/) {
      $stophere=1;
    }
  }

  # OK now we have the optimal lambda. Interpolate the intersections sliced by
  # the plane at the optimal height
  my $sx=0, $sy=0;
  for $i ($#xhis) {
    $sx += $lam*$xhis[$i] + $mal*$xlos[$i];
    $sy += $lam*$yhis[$i] + $mal*$ylos[$i];
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

  my $inside_am = ($dh < $area/$n/$MAGIC_NUMBER);
  my $covar = Math::MatrixReal->new_from_rows( [ [$varx, $cvar],
						 [$cvar, $vary] ]);
  my $icov = $covar->inverse();
  my $ellipse_stat = ($icov->element(1,1) * $avgx*$avgx
                    + $icov->element(2,2) * $avgy*$avgy
                  + 2*$icov->element(1,2) * $avgx*$avgy);
  my $inside_e39 = ($ellipse_stat < 1.000);
  my $inside_e90 = ($ellipse_stat < 4.605);
  my $inside_e95 = ($ellipse_stat < 5.991);

  # if ($n >= 50) {
  # $total_count++;
  # $total_inside += $inside_e90;
  # $total_avg = $total_inside / $total_count;
  # print STDERR "INSIDE_AM\t$n\t$total_inside\t$total_count\t$total_avg\n";
  # }


  return (join ',', 'POLY', $n, $area, $avgx, $avgy, $dh, $z,
	            $lam, $varx, $vary, $cvar, 
	            $icov->element(1,1), $icov->element(2,2), $icov->element(1,2),
	            $inside_am, $ellipse_stat, $inside_e39, $inside_e90, $inside_e95, $ambig, "\n");
}




sub center_spread {
  # compute the area of a 1-sigma error ellipse
  my ($n,$sx,$sy,$sxx,$syy,$sxy)=(0,0,0,0,0,0);
  while (@_) {
    my $x = shift;
    my $y = shift;
    $n++;
    $sx += $x;
    $sy += $y;
    $sxx += $x*$x;
    $syy += $y*$y;
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

sub hourglass {
  my @axyzs;
  for $planei (0..$nplanes) {
    $phi = $planei/$nplanes;
    $plo = 1 - $phi;

    @xys = ();
    for $iid (@_) {
      $h = $hsh->{$iid};
      $x = $phi * $h->{xhi} + $plo * $h->{xlo};
      $y = $phi * $h->{yhi} + $plo * $h->{ylo};
      $z = $phi * $h->{zhi} + $plo * $h->{zlo};
      push @xys, $x, $y;
    }
    ($avgx, $avgy, $area) = center_spread(@xys);
    $zstr = sprintf "%.1f", $z;
    push @axyzs, (join ',', $area, $avgx, $avgy, $zstr, $phi);
  }
  my @saxyzs = sort {$a<=>$b} @axyzs;
  return $saxyzs[0];
}


sub random_images {
  my $n = shift;
  my %unused = ();
  my %used = ();
  for $iid (@iids) { $unused{$iid} = 1 }
  while (keys %used < $n) {
    my $nleft = keys %unused;
    my $i = rand()*$nleft;
    my $iid = (keys %unused)[$i];
    delete $unused{$iid}; # so we don't sample it again
    $used{$iid}++;
  }
  return keys %used;
}


#########################################################

select STDOUT; $| = 1; # autoflush
$pi = atan2(1,1)*4;
$nplanes = 200; # plus 1!
for $i (0..$nplanes) { push @lams, $i/$nplanes }

# 36N 117.5W 1700 in UTM11N
$truthx =  454936.0104;
$truthy = 3984064.0306;
$truthz = 1700;

$MAGIC_NUMBER = 9.8;

@iids = ();
while (<>) {
  chomp;
  @ary = split /,/;
  $iid = shift @ary;
  push @iids, $iid;
  for $key (qw(xhi yhi zhi  xmd ymd zmd  xlo ylo zlo)) {
    $hsh->{$iid}->{$key} = shift @ary;
    if ($key =~ /x/)    { $hsh->{$iid}->{$key} -= $truthx }
    if ($key =~ /y/)    { $hsh->{$iid}->{$key} -= $truthy }
    #if ($key =~ /[xy]/) { $hsh->{$iid}->{$key} *= 100000 }
    #if ($key eq 'zhi')  { $hsh->{$iid}->{$key}  =   10.0 }
    #if ($key eq 'zlo')  { $hsh->{$iid}->{$key}  =  -10.0 }
    if ($key =~ /z/)    { $hsh->{$iid}->{$key} -= $truthz }
  }
}


print join ',', qw(ALG N AREA AVGX AVGY DH Z
                   L_MIN VARX VARY CVARXY INVX INVY INVXY
                   INA ELL INE39 INE90 INE95 AMBIG), "\n";
@ns = ();
for ($n=4;   $n<100;  $n++)  { push @ns, $n } # every $n=4..99
for ($n=100; $n<1000; $n+=5) { push @ns, $n } # 100...995 by 5s
for $n (@ns) {
  for (1..100) {
    @sample = random_images($n);
#    $hg  = hourglass(@sample);
#    $hgb = hourglass_brute(@sample);
    $hgp = hourglass_poly(@sample);
#    $poly_gain = $hgb - $hgp;;
#    $hgp .= ",$poly_gain";
#    print 'OLDE,', $n, ',', $hg,  "\n";
#    print 'BRUT,', $n, ',', $hgb, "\n";
    print $hgp;
  }
}

# and finally 1000, only one sample possible
$n = 1000;
#print 'OLDE,', $n, ',', hourglass(@iids),       "\n";
#print 'BRUT,', $n, ',', hourglass_brute(@iids), "\n";
print hourglass_poly(@iids);
