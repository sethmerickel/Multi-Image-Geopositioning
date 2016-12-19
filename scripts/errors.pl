#! /usr/bin/perl


sub sample {
  my $dx = shift;
  my $dy = shift;
  my $dz = shift;
  $s{x} += $dx;
  $s{y} += $dy;
  $s{z} += $dz;
  $s{xx} += $dx*$dx;
  $s{yy} += $dy*$dy;
  $s{zz} += $dz*$dz;
  $s{xy} += $dx*$dy;
  $s{xz} += $dx*$dz;
  $s{yz} += $dy*$dz;
  $nnn++;
}

sub clear {
  %s = ();
  $nnn = 0;
}

sub covar {
  my $l1 = shift;
  my $l2 = shift;
  return ($s{$l1.$l2} - $s{$l1}*$s{$l2}/$nnn)/($nnn-1);
}

sub absavg {
  my $sum=0;
  my $n=@_;
  for my $x (@_) { $sum += abs($x) }
  return $sum/$n;
}


print (join ',', qw(N AVGDH AVGDV 
                           DH39 DH68 DH95 DH99 RH95 RH99
                           DZ68 DZ95 DZ99 RZ95 RZ99
                           AREA A_N_DH39 A_N_DH68
                           E39 E90 E95
), "\n");

$K = 100;

while (<>) {
  #($n,$a,$dx,$dy,$dh,$dz,$e,$ine39,$ine90,$ine95) = (split /,/)[1,2,3,4,5,6,15,16,17,18];
  ($n,$x,$y,$z) = (split /,/)[1,2,3,4];
  $dx = $x-455000;
  $dy = $y-3984000;
  $dz = $z-1700;
  $dh = sqrt($dx*$dx + $dy*$dy);
  next unless $n > 0;
  push @dxs, $dx;
  push @dys, $dy;
  push @dhs, $dh;
  push @dzs, $dz;
  push @es,  $e;
  $suma += $a;
  sample($dx,$dy,$dz);
  $cte39 += $ine39;
  $cte90 += $ine90;
  $cte95 += $ine95;
  if ($nnn == $K) {
    # report and clear out
    @sdxs = sort {$a<=>$b} @dxs;
    @sdys = sort {$a<=>$b} @dys;
    @sdhs = sort {$a<=>$b} @dhs;
    @sdzs = sort {$a<=>$b} @dzs;
    @ses  = sort {$a<=>$b} @es;
    $dh39 = $sdhs[.38*$K]; 
    $dh68 = $sdhs[.67*$K]; $dz68 = $sdzs[.67*$K];
    $dh95 = $sdhs[.94*$K]; $dz95 = $sdzs[.94*$K];
    $dh99 = 0.7*$sdhs[.98*$K] + 0.3*$sdhs[.99*$K];
    $dz99 = 0.7*$sdzs[.98*$K] + 0.3*$sdzs[.99*$K];
    $rh95 = $dh95 / $dh68; $rz95 = $dz95 / $dz68;
    $rh99 = $dh99 / $dh68; $rz99 = $dz99 / $dz68;

    $avgdh = absavg(@dhs);
    $avgdv = absavg(@dzs);
    $avga = $suma / $nnn;
    print (join ',', $n, $avgdh, $avgdv,
                                $dh39, $dh68, $dh95, $dh99, $rh95, $rh99,
                                $dz68, $dz95, $dz99, $rz95, $rz99,
                                $avga, $avga/$n/$dh39, $avga/$n/$dh68,
                                $cte39/$K, $cte90/$K, $cte95/$K,
	   "\n");
    @dxs = @dys = @dhs = @dzs = ();
    $suma = 0;
    $cte39 = $cte90 = $cte95 = 0;
    clear();
  }
}
