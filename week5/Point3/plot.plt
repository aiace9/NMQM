set output "onda_s.ps"
set terminal postscript color
set size 1.,1.
set title 'onda s'
plot[0:10][] 's-wfc-3.out' u 1 : (-$3) w l, 's-wfc-4.out' u 1 : ($3) w l

set output "onda_p.ps"
set terminal postscript color
set size 1.,1.
set title 'onda p'
plot[0:10][] 'p-wfc-3.out' u 1 : (-$3) w l, 'p-wfc-4.out' u 1 : ($3) w l