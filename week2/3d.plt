set output "data_p.ps"
set terminal postscript portrait
set size 1.,1.
set yrange[-0.5:0.5]
set xrange[0:30]
set xlabel "radius"
set grid
show grid
plot 'pot.out' u 1 : 2 w l, 'wfc.out' u 1 : ($3)**2 w l, 'wfc.out' u 1 : 2 w l

