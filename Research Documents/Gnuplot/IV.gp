#!/usr/bin/env gnuplot --persist
set terminal postscript eps enhanced color "Helvetica" 24
#set terminal pdf enhanced fname 'Helvetica' fsize 10
#set terminal aqua
#set datafile separator ","

set style line 1 lw 4
set style line 2 lw 4
set style line 3 lw 4
set style line 4 lw 4

# Pattern Radious (in cm)
r1 = 0.05
r2 = 0.1
# Thus, pattern area becomes,
area1 = r1**2 * 3.14
area2 = r2**2 * 3.14

# Thickness (in cm)
thick = 300e-7 

set output 'JE.eps'
set log y
set xlabel "Field (MV/cm)"
set ylabel "Current Density (A/cm^{2})"
set format y "10^{%L}"
VI_TITLE = sprintf("Dark Current Density of a-SiN_{x}:H (thickness: %d nm)", thick*1e7)
set title VI_TITLE
set key bottom right
show margin 

# Normalizing Data to E/J format while plotting
plot "1mm.csv" using ($1/thick/1e6):($2/(area1)) title "1mm Pos. 1" w l ls 1, \
    "1mm.csv" using ($1/thick/1e6):($4/(area1)) title "1mm Pos. 2" w l ls 2, \
    "2mm.csv" using ($1/thick/1e6):($2/(area2)) title "2mm Pos. 1" w l ls 3, \
    "2mm.csv" using ($1/thick/1e6):($4/(area2)) title "2mm Pos. 2" w l ls 4


set output 'IV.eps'
set log y
set xlabel "Bias (V)"
set ylabel "Current (A)"
set format y "10^{%L}"
VI_TITLE = sprintf("Dark Current of a-SiN_{x}:H (thickness: %d nm)", thick*1e7)
set title VI_TITLE
set key bottom right
show margin 

# Normalizing Data to E/J format while plotting
plot "1mm.csv" using ($1):($2) title "1mm Pos. 1" w l ls 1, \
    "1mm.csv" using ($1):($4) title "1mm Pos. 2" w l ls 2, \
    "2mm.csv" using ($1):($2) title "2mm Pos. 1" w l ls 3, \
    "2mm.csv" using ($1):($4) title "2mm Pos. 2" w l ls 4


