#set output "stress_plot.eps"
#set terminal postscript eps enhanced color size 3.5,2.5

set grid
set xlabel "Applied stress"
set ylabel "Kink force [GPa]"
set yrange [-4:4]
set xrange [0:1]

a = 3.18						# Angstr
b = a*sqrt(3)/2					# Angstr
mu = 160# e9# * 6.241506e-12		# Shear modulus in [eV/A^3] units:
h = a*sqrt(6)/3					# Angstr
poissonRatio = 0.28
poissonFactor = (1 + poissonRatio) / (1 - poissonRatio)
w0_1 = 4.1377
w0_2 = 3.95369
wp = 3.31418
peierls_stress = 3.2

w1(sigma) = w0_1*log(1/sigma)
w2(sigma) = w0_2*log(1/sigma)+wp

f(ksp) = mu*h*h/(8*pi*ksp*ksp)*poissonFactor

kink_force1(sigma) = f(w1(sigma)*b)-sigma*peierls_stress
kink_force2(sigma) = f(w2(sigma)*b)-sigma*peierls_stress

plot kink_force1(x), kink_force2(x)
