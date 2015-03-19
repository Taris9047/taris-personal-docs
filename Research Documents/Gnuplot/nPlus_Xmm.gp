#!/usr/bin/env gnuplot --persist
set terminal postscript eps enhanced "Helvetica" 24
#set terminalinal pdf enhanced fname 'Helvetica' fsize 14
#set terminal aqua
set datafile separator ","
set style line 1 lw 4

#
# Setting up shadow mask dimension data.
# W, L in cm
W = 2
L = 0.1 # Usually, this needs to be changed!
# nPlus Thickness (cm).
thick = 50e-7

# Output Filename
OUT_FILENAME = 'nPlus1mm.eps'

# Plot Title
PLOT_TITLE = sprintf("nPlus from Plasmatherm, Thickness = %d nm", thick*1e7)

set output OUT_FILENAME
#set log y
set xlabel "Voltage Bias (V)"
set ylabel "Current (A)"
set format y "%1.1E"
#set log y
#set format y "10^{%L}"
#set yrange [-2e-6:2e-6]
set title PLOT_TITLE
set key bottom right
show margin 
#
# Calculating Resistivity/Conductivity
f(x) = a*x
fit f(x) "1mm.csv" using ($1):(($2+$3+$4+$5)/4) via a 
# Now, A contains the Conductance. So, the resistance is,
Resistance = 1/a
# Therefore, Conductance is ...
Conduct = L/(Resistance*W*thick)

set label "{/Symbol s} = %.3f S/cm", Conduct at graph 0.05, graph 0.95
plot "1mm.csv" using ($1):($2) title "Position 1" w l, \
	 "1mm.csv" using ($1):($3) title "Position 2" w l, \
	 "1mm.csv" using ($1):($4) title "Position 3" w l, \
	 "1mm.csv" using ($1):($5) title "Position 4" w l

	
