# Gnuplot MS-Windows 32 bit version 4.0
# Plots the energy distribution of primary particles simulated by the
# example main programs 'penslab.f', 'pencyl.f' and 'penmain.f'

# unset mouse

set xzeroaxis
set yrange [:]

set title 'Energy distribution of primary particles'
set xlabel "E (eV)"
set ylabel "p(E) (1/eV)"
plot 'psource.dat' u 1:2 t 'source distribution' w l 1, \
     'psource.dat' u 1:3 t 'MC spectrum - low  ' w l 2, \
     'psource.dat' u 1:4 t 'MC spectrum - high ' w l 3
pause -1
