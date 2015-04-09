set output "data_p.ps"
set terminal postscript portrait
set size 1.,1.
plot[0:10][] 's-wfc.out' u 1 : (-$3*sqrt(4*pi)) w l, 'wfc.out' u 1:3 w l
