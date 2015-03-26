set output "data_p.ps"
set terminal postscript portrait
set size 1.,1.
set xlabel "Energy"
set ylabel "total scattering section"
set grid
show grid
plot 'data.out' u 1 : 2 w p
