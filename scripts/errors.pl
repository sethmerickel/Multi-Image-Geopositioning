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


print (join ',', qw(N DH39 DH68 DH95 DH99 RH95 RH99
                           DZ68 DZ95 DZ99 RZ95 RZ99
                           AREA A_N_DH39 A_N_DH68), "\n");

while (<>) {
  ($n,$a,$dx,$dy,$dh,$dz) = (split /,/)[1,2,3,4,5,6];
  next unless $n > 0;
  push @dxs, $dx;
  push @dys, $dy;
  push @dhs, $dh;
  push @dzs, $dz;
  $suma += $a;
  sample($dx,$dy,$dz);
  if ($nnn == 100) {
    # report and clear out
    @sdxs = sort {$a<=>$b} @dxs;
    @sdys = sort {$a<=>$b} @dys;
    @sdhs = sort {$a<=>$b} @dhs;
    @sdzs = sort {$a<=>$b} @dzs;
    $dh39 = $sdhs[38]; 
    $dh68 = $sdhs[67]; $dz68 = $sdzs[67];
    $dh95 = $sdhs[94]; $dz95 = $sdzs[94];
    $dh99 = 0.7*$sdhs[98] + 0.3*$sdhs[99];
    $dz99 = 0.7*$sdzs[98] + 0.3*$sdzs[99];
    $rh95 = $dh95 / $dh68; $rz95 = $dz95 / $dz68;
    $rh99 = $dh99 / $dh68; $rz99 = $dz99 / $dz68;

    $varx = covar('x','x');
    $vary = covar('y','y');
    $varz = covar('z','z');
    $cvxy = covar('x','y');
    $cvxz = covar('x','z');
    $cvyz = covar('y','z');

    $avga = $suma / $nnn;
    print (join ',', $n, $dh39, $dh68, $dh95, $dh99, $rh95, $rh99,
                                $dz68, $dz95, $dz99, $rz95, $rz99,
                                $avga, $avga/$n/$dh39, $avga/$n/$dh68,"\n");
    @dxs = @dys = @dhs = @dzs = ();
    $suma = 0;
    clear();
  }
}
