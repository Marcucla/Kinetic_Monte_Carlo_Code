#!/opt/local/bin/gnuplot

set terminal pngcairo size 480,362 enhanced font 'Verdana,10'
set output 'kink_pair_separation.png'

set yrange [0:20]
set xrange [0:1]
set xlabel "Resolved shear stress / Peierls stress"
set ylabel "Kink pair separation (b)"
set grid

kpsep_fit(sigma)=w0*log(1/sigma)+wp
fit kpsep_fit(x) 'Fe_data.txt' using ($1/1e3):($2) via w0,wp

plot 'Fe_data.txt' using ($1/1e3):($2) w linespoints, kpsep_fit(x) title "w_0 * log(1/sigma) + w_p"