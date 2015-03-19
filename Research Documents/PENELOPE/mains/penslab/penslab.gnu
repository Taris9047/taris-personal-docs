# Gnuplot MS-Windows 32 bit version 4.0
# Plots results from 'penslab.f'

# unset mouse

set xzeroaxis
set yrange [:]

set title 'Depth dose distribution'
set xlabel "depth (cm)"
set ylabel "dose (eV/(g/cm**2))"
plot 'psddose.dat' u 1:2:3 w errorbars 2, \
     'psddose.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Depth distribution of deposited charge'
set xlabel "depth (cm)"
set ylabel "charge density (e/cm)"
plot 'psdchar.dat' u 1:2:3 w errorbars 2, \
     'psdchar.dat' u 1:2  notitle w histeps 1
pause -1

set title 'Path length distribution of transmitted primary particles'
set xlabel "path length (cm)"
set ylabel "PDF (1/cm)"
plot 'pstltr.dat' u 1:2:3 w errorbars 2,\
     'pstltr.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Path length distribution of backscattered primary particles'
set xlabel "path length (cm)"
set ylabel "PDF (1/cm)"
plot 'pstlbk.dat' u 1:2:3 w errorbars 2,\
     'pstlbk.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Path length distribution of absorbed primary particles'
set xlabel "path length (cm)"
set ylabel "PDF (1/cm)"
plot 'pstlab.dat' u 1:2:3 w errorbars 2,\
     'pstlab.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Angular distribution of emerging electrons'
set xlabel "angle (deg)"
set ylabel "PDF (1/sr)"
plot 'psanel.dat' u 1:2:3 w errorbars 2,\
     'psanel.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Angular distribution of emerging photons'
set xlabel "angle (deg)"
set ylabel "PDF (1/sr)"
plot 'psanga.dat' u 1:2:3 w errorbars 2,\
     'psanga.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Angular distribution of emerging positrons'
set xlabel "angle (deg)"
set ylabel "PDF (1/sr)"
plot 'psanpo.dat' u 1:2:3 w errorbars 2,\
     'psanpo.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of transmitted electrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pseneltr.dat' u 1:2:3 w errorbars 2,\
     'pseneltr.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of backscattered electrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'psenelbk.dat' u 1:2:3 w errorbars 2,\
     'psenelbk.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of transmitted photons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'psengatr.dat' u 1:2:3 w errorbars 2,\
     'psengatr.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of backscattered photons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'psengabk.dat' u 1:2:3 w errorbars 2,\
     'psengabk.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of transmitted positrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'psenpotr.dat' u 1:2:3 w errorbars 2,\
     'psenpotr.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of backscattered positrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'psenpobk.dat' u 1:2:3 w errorbars 2,\
     'psenpobk.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Deposited energy distribution'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'psedepm.dat' u 1:2:3 w errorbars 2, \
     'psedepm.dat' u 1:2 notitle w histeps 1
pause -1

