#!/usr/bin/env gnuplot --persist
set terminal postscript eps enhanced color "Helvetica" 24
#set terminal pdf enhanced fname 'Helvetica' fsize 10
#set terminal aqua
set style line 1 lw 4

set datafile separator ","
set output 'SiNx_CV.eps'


# Pattern Radious (in cm)
r1 = 0.05
r2 = 0.1
# Thus, pattern area becomes,
area1 = r1**2 * 3.14
area2 = r2**2 * 3.14

# Thickness (in cm)
thick = 190e-7

set xlabel "Voltage (V)"
set ylabel "Normalized Capacitance (A. U.)"
#set format y "%1.0E"
set title "C-V Sweep of a-SiN_{x}"
set key top right
show margin


plot "1mm1MHz_40.csv" using ($1):($4) title "1 MHz Position 1" w l, \
    "2mm1MHz_40.csv" using ($1):($4) title "1 MHz Position 2" w l

	
