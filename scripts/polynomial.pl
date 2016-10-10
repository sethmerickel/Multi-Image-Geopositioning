#! /usr/bin/perl
use FindBin;
use lib $FindBin::Bin;
use Hourglass ':all';


select STDOUT; $| = 1; # autoflush
#$nplanes = 200; # plus 1!

# 36N 117.5W 1700 in UTM11N
$truthx =  454936.0104;
$truthy = 3984064.0306;
$truthz = 1700;


open HIMIDLO, "$ARGV[0]";
@lines = (<HIMIDLO>);
$hsh = parse_himidlo($truthx, $truthy, $truthz, @lines);
@iids = keys %$hsh;

print join ',', qw(ALG N AREA AVGX AVGY DH Z
                   L_MIN VARX VARY CVARXY INVX INVY INVXY
                   ELL INE39 INE90 INE95 AMBIG), "\n";
@ns = ();
for ($n=4;   $n<100;  $n++)  { push @ns, $n } # every $n=4..99
for ($n=100; $n<1000; $n+=5) { push @ns, $n } # 100...995 by 5s
for $n (@ns) {
  for (1..100) {
    @sample = random_images($n, @iids);
#    $hg  = hourglass(@sample);
#    $hgb = hourglass_brute(@sample);
    $hgp = hourglass_poly($hsh, @sample);
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
print hourglass_poly($hsh, @iids);
