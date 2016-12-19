#! /usr/bin/perl

for $f (@ARGV) { # output files from migpoly.pl
  open MIGPOLY, $f;
  while (<MIGPOLY>) {
    next if /ALG,N/; # skip the header line
    ($n, $mx, $my, $mz, $px, $py, $pz) = (split /,/)[1, 2,3,4, 15,16,18];
    $mx -= 455000; $my -= 3984000; $mz -= 1700;
    $px -= 455000; $py -= 3984000; $pz -= 1700;
    $mh = sqrt($mx*$mx + $my*$my); $mv = abs($mz);
    $ph = sqrt($px*$px + $py*$py); $pv = abs($pz);
    $stat->{$f}->{$n}->{smh} += $mh; $stat->{$f}->{$n}->{smv} += $mv;
    $stat->{$f}->{$n}->{sph} += $ph; $stat->{$f}->{$n}->{spv} += $pv;
    $stat->{$f}->{$n}->{num} += 1;
  }
  close MIGPOLY;
}

@hary1 = ('');
for $f (@ARGV) { push @hary1, $f; push @hary1, ('')x3; }
push @hary1, '', '';
for $f (@ARGV) { push @hary1, $f; push @hary1, ('')x3; }
print (join ',', @hary1, "\n");

@hary2 = ('N');
for $f (@ARGV) { push @hary2, qw(MDH PDH MDV PDV) }
push @hary2, '', 'N';
for $f (@ARGV) { push @hary2, qw(MDH PDH MDV PDV) }
print (join ',', @hary2, "\n");

@ns = sort {$a<=>$b} keys %{$stat->{$ARGV[0]}};
for $n (@ns) {
  @pary = ($n);
  for $f (@ARGV) {
    $s = $stat->{$f}->{$n};
    $fact = 1.0/$s->{num};
    push @pary, $s->{smh} * $fact;
    push @pary, $s->{sph} * $fact;
    push @pary, $s->{smv} * $fact;
    push @pary, $s->{spv} * $fact;
  }
  push @pary, '', $n;
  for $f (@ARGV) {
    $s = $stat->{$f}->{$n};
    $fact = 1.0/$s->{num} * sqrt($n);
    push @pary, $s->{smh} * $fact;
    push @pary, $s->{sph} * $fact;
    push @pary, $s->{smv} * $fact;
    push @pary, $s->{spv} * $fact;
  }
  print (join ',', @pary, "\n");

  last if $n >= 200;
}
