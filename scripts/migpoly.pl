#! /usr/bin/perl
use FindBin;
use lib $FindBin::Bin;
use Hourglass ':all';

use Getopt::Std;
%opt = (n=>100, m=>0);
getopts('n:m:hqio', \%opt);


select STDOUT; $| = 1; # autoflush
#$nplanes = 200; # plus 1!

# near 36N 117.5W 1700 in UTM11N
#$truthx =  455000;
#$truthy = 3984000;
#$truthz =    1700;
$truthx = $truthy = $truthz = 0;

@lines = (<>); # slurp up projector log
$hsh = parse_projector($truthx, $truthy, $truthz, @lines);

# if half or quarter is selected, delete southern/western looks
if ($opt{h} || $opt{q}) { north_images($hsh) }
if            ($opt{q}) {  east_images($hsh) }
if ($opt{i})            { inner_images($hsh) }
if ($opt{o})            { outer_images($hsh) }
@iids = sort keys %$hsh;
$nimg = @iids;

# simulate measurement error (default 0)
add_meas_err($hsh, $opt{m});

print join ',', qw(
ALG N X Y Z REF VARX CVARXY CVARXZ VARY CVARYZ VARZ
ALG N AREA AVGX AVGY DH Z
L_MIN VARX VARY CVARXY INVX INVY INVXY
ELL INE39 INE90 INE95 AMBIG), "\n";
@ns = ();
for ($n=4;   $n<100  && $n<$nimg; $n+=1) { push @ns, $n } # every $n=4..99
for ($n=100; $n<1000 && $n<$nimg; $n+=5) { push @ns, $n } # 100...995 by 5s
for $n (@ns) {
  for (1..$opt{n}) {
    @sample = random_images($n, @iids);

    $gp0 = mkmat(3,1, 1,1,1);
    ($gp, $cov, $refvar) = wvmig($hsh, $gp0, @sample);
    $toprint = join ',', 'MIG', $n, flatten($gp), $refvar,
                           (flatten($cov))[0,1,2,4,5,8], '';

    $toprint .= hourglass_poly($hsh, @sample);

    print $toprint;
    #print "\n";
  }
}

# and finally all the images, only one sample possible
$n = $nimgs;
$gp0 = mkmat(3,1, 1,1,1);
($gp, $cov, $refvar) = wvmig($hsh, $gp0, @iids);
$toprint = join ',', 'MIG', $n, flatten($gp), $refvar,
                    (flatten($cov))[0,1,2,4,5,8], '';
$toprint .= hourglass_poly($hsh, @iids);

print $toprint;
#print "\n";
