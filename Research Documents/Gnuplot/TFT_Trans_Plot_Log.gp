#!/usr/bin/env gnuplot --persist
set terminal postscript eps enhanced "Helvetica" 24
#set terminal pdf enhanced fname 'Helvetica' fsize 14
#set terminal aqua
#set datafile separator ","
set style line 1 lw 4


#
# Setting up TFT dimension data.
# W/L in cm
W = 200e-4
L = 15e-4 # Usually, this needs to be changed!
#
# Output Filename
OUT_FILENAME = 'SampleTrans.eps'
# Input Filename
IN_FILENAME = 'SampleTrans.csv'

# Plot Title
PLOT_TITLE = sprintf("Sample a-Si:H TFT, W/L = %d/%d", W*1e4, L*1e4)

set output OUT_FILENAME
set log y
set xlabel "Gate Bias (V)"
set ylabel "Drain Current (A)"
set format y "10^{%L}"
set title PLOT_TITLE
set key bottom right
show margin 
plot IN_FILENAME using 1:2 title "V_{ds} = 0.1 V" w l, \
	IN_FILENAME using 1:3 title "V_{ds} = 1.0 V" w l, \
	IN_FILENAME using 1:4 title "V_{ds} = 10 V" w l, \
    IN_FILENAME using 1:5 title "V_{ds} = 20 V" w l, \
    IN_FILENAME using 1:6 title "V_{ds} = 30 V" w l
