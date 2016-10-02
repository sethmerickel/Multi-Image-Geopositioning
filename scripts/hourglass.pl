#! /usr/bin/perl
use Math::MatrixReal;

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
  if ($area < 0) {
    $stophere=1;
  }
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
    push @axyzs, (join ',', $area, $avgx, $avgy, $z);
  }
  my @saxyzs = sort {$a<=>$b} @axyzs;
  return $saxyzs[0];
}


sub random_images {
  my $n = shift;
  my %u = ();
  while (keys %u < $n) {
    my $i = rand()*@iids;
    $u{$iids[$i]}++;
  }
  return keys %u;
}


#########################################################

$pi = atan2(1,1)*4;
$nplanes = 200; # plus 1!

@iids = ();
while (<>) {
  chomp;
  @ary = split /,/;
  $iid = shift @ary;
  push @iids, $iid;
  for $key (qw(xhi yhi zhi  xmd ymd zmd  xlo ylo zlo)) {
    $hsh->{$iid}->{$key} = shift @ary;
    if ($key =~ /x/)    { $hsh->{$iid}->{$key} -= -117.5 }
    if ($key =~ /y/)    { $hsh->{$iid}->{$key} -=   36.0 }
    if ($key =~ /[xy]/) { $hsh->{$iid}->{$key} *= 100000 }
    if ($key eq 'zhi')  { $hsh->{$iid}->{$key}  =   10.0 }
    if ($key eq 'zlo')  { $hsh->{$iid}->{$key}  =  -10.0 }
  }
}


for $n (5..$#iids) {
  @sample = random_images($n);
  print $n, ',', hourglass(@sample), "\n";
}
