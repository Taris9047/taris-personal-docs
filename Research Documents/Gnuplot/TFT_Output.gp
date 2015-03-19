#!/usr/bin/env gnuplot --persist
set terminal postscript eps enhanced "Helvetica" 24
#set terminal pdf enhanced fname 'Helvetica' fsize 14
#set terminal aqua
set datafile separator ","
set style line 1 lw 4


clear
#
# Setting up TFT dimension data.
# W/L in cm
W = 200e-4
L = 15e-4 # Usually, this needs to be changed!
#
# Output Filename
OUT_FILENAME = 'SampleOutput.eps'
# Input Filename
IN_FILENAME = 'SampleOutput.csv'

# Plot Title
PLOT_TITLE = sprintf("Sample a-Si:H TFT, W/L = %d/%d", W*1e4, L*1e4)

set output OUT_FILENAME
set xlabel "Drain Bias (V)"
set ylabel "Drain Current (A)"
set format y "%1.1E"
set title PLOT_TITLE
set key top left
show margin 
plot IN_FILENAME using 1:2 title "V_{gs} = 10 V" w l, \
	IN_FILENAME using 1:3 title "V_{gs} = 15 V" w l, \
	IN_FILENAME using 1:4 title "V_{gs} = 20 V" w l, \
    IN_FILENAME using 1:5 title "V_{gs} = 25 V" w l, \
    IN_FILENAME using 1:6 title "V_{gs} = 30 V" w l
