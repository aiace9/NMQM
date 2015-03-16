set output "data_p.ps"
set terminal postscript portrait
set size 1.,1.
set yrange[-0.5:1.3]
set xrange[0:5]
set xlabel "radius"
set ylabel "prob"
set grid
show grid
Z=1
a0= 1
spi=sqrt(pi)
plot 'wfc.out' u 1 : ($3)**2 w l, 'wfc.out' u 1 : 2 w l, ((Z/a0)**(3/2))*(1/spi)*exp(-Z*x/a0)
