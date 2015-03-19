#!/usr/bin/env gnuplot --persist
set terminal postscript eps enhanced "Helvetica" 24
#set terminalinal pdf enhanced fname 'Helvetica' fsize 14
#set terminal aqua
set datafile separator ","
set style line 1 lw 4


#
# Setting up TFT dimension data.
# W/L in cm
W = 200e-4
L = 15e-4 # Usually, this needs to be changed!
# Dielectric Thickness (cm).
tox = 300e-7
# Dielectric Constant of insulator.
eps_ox = 6.3
#
# Input/Output Filename
IN_FILENAME = 'SampleTrans.csv'
OUT_FILENAME = 'SampleFitLin.eps'

# Plot Title
PLOT_TITLE = sprintf("Sample a-Si:H TFT, W/L = %d/%d", W*1e4, L*1e4)

set output 
#set log y
set xlabel "Gate Bias (V)"
set ylabel "Drain Current (A)"
set yrange [-1E-12:]
set format y "%1.1E"
set title PLOT_TITLE
set key bottom right
show margin 
# defining a linear fit function for curve fitting.
f(x) = a*x + b
# Fitting the linear portion of source graph. The fitting range requires adjustment everytime.
# Change [10:20] part according to your TFT characteristic.
fit [10:20] f(x) IN_FILENAME using ($1):($3) via a, b

#
# Calculating Field effect mobility
# Calculating Cox (or, Csinx)
Cox = eps_ox*8.854e-14/tox
# Finally, Mu is,
Mu = a*(L/(W*Cox))
#

set label "{/Symbol m}_{EFF} = %1.2f cm^{2}/V{/Symbol \327}s", Mu at graph 0.05, graph 0.95
set label "V_{Th} = %1.2f V", -b/a at graph 0.05, graph 0.85
plot IN_FILENAME using ($1):($3) title "V_{ds} = 1 V" w l, \
	f(x) title "Linear Fit" w l
