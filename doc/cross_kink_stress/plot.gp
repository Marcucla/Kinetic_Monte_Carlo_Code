set output "stress_plot.eps"
set terminal postscript eps enhanced color size 3.5,2.5

set grid
set xlabel "Kink separation [b]"
set ylabel "Shear stress [GPa]"
#set yrange [-3:3]

a = 3.18						# Angstr
b = a*sqrt(3)/2					# Angstr
mu = 160# e9# * 6.241506e-12		# Shear modulus in [eV/A^3] units:
h = a*sqrt(6)/3					# Angstr
poissonRatio = 0.28
poissonFactor = (1 + poissonRatio) / (1 - poissonRatio)

plot \
	'cross_kink_60_stress_0.5b.dat' using 1:($2/1e9) title "integral" w l lw 2, \
	'cross_kink_60_stress_0.5b_old.dat' using 1:($2/1e9) title "center" w l lw 2, \
	'cross_kink_60_stress_0.5b_0.dat' using 1:($2/1e9) title "corner1" w l lw 2, \
	'cross_kink_60_stress_0.5b_1.dat' using 1:($2/1e9) title "corner2" w l lw 2
	