#!/opt/local/bin/gnuplot

set terminal pngcairo size 480,362 enhanced font 'Verdana,10'
set output 'kink_pair_energy.png'

set yrange [-3:+3]
set xrange [0:1]
set xlabel "{/Symbol s}/{/Symbol s}_p"
set ylabel "Energy (eV)"
set grid

kT = 300 * 8.6173324e-5			# Boltzmann constant times temperature (in eV units).
a = 3.18						# Angstr
b = a*sqrt(3)/2					# Angstr
mu = 160e9 * 6.241506e-12		# Shear modulus in [eV/A^3] units:
h = a*sqrt(6)/3					# Angstr
poissonRatio = 0.28
peierls_stress = 3.2e9			# Pa
deltaH0 = 1.75					# eV
w0 = 4.1377						# b

poissonFactor = (1 + poissonRatio) / (1 - poissonRatio)

deltaH(sigma) = deltaH0 * (1.0 - sqrt(sigma))**1.25
w(sigma) = w0 * log(1 / sigma)
elasticInteraction(sigma) = -mu * b * h * h / (8 * pi * w(sigma)) * poissonFactor
kp_energy(sigma) = deltaH(sigma) + elasticInteraction(sigma)

plot deltaH(x) title "{/Symbol D}H={/Symbol D}H_0 (1-sqrt({/Symbol s}/{/Symbol s}_p))^{1.25}", elasticInteraction(x) title "Elastic interaction", kp_energy(x) title "Kink-pair activation energy"