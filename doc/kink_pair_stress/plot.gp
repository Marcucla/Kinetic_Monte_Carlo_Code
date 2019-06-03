set output "stress_plot.eps"
set terminal postscript eps enhanced color size 3.5,2.5

set grid
set xlabel "Kink separation [b]"
set ylabel "Shear stress [GPa]"
set yrange [-3:3]

a = 3.18						# Angstr
b = a*sqrt(3)/2					# Angstr
mu = 160# e9# * 6.241506e-12		# Shear modulus in [eV/A^3] units:
h = a*sqrt(6)/3					# Angstr
poissonRatio = 0.28
poissonFactor = (1 + poissonRatio) / (1 - poissonRatio)

f(ksp) = mu*h*h/(8*pi*ksp*ksp)*poissonFactor/b/h

plot \
	'kink_pair_stress_sharp_0.5b.dat' using 1:($2/1e9) title "Sharp kinks" w l lw 2, \
	'kink_pair_stress_beveled_0.5b.dat' using 1:($2/1e9) title "Beveled kinks" w l lw 2, \
	f(x) title "Hirth & Lothe force formula"
