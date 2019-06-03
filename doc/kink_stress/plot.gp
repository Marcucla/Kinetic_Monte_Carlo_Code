set output "stress_plot.eps"
set terminal postscript eps enhanced color size 3.5,2.5

set grid
set xlabel "Position along screw direction [b]"
set ylabel "Shear stress [GPa]"
plot 'sharp_kink_0.5b.dat' using 1:($2/1e9) title "Sharp kink (case 1)" w l lw 2, \
	'wide_kink_0.5b.dat' using 1:($2/1e9) title "Beveled kink (case 2)"  w l lw 2 lc 3, \
	'wide_kink_on_0.5b.dat' using 1:($2/1e9) title "Beveled kink (case 3)"  w l lw 2 lc 4


set output "stress_plot_case3.eps"
set size 1,0.6
plot 'wide_kink_on_0.5b.dat' using 1:($2/1e9) title "Beveled kink (case 3)"  w l lw 2 lc 4
	