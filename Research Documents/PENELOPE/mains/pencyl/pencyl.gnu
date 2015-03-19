# Gnuplot MS-Windows 32 bit version 4.0
# Plots results from 'pencyl.f'

# unset mouse
set zero 1.0e-60

set xzeroaxis
set yrange [:]

set title 'Depth dose distribution'
set xlabel "depth (cm)"
set ylabel "dose (eV/(g/cm**2))"
plot 'pcddose.dat' u 1:2:3 w errorbars 2, \
     'pcddose.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Depth distribution of deposited charge'
set xlabel "depth (cm)"
set ylabel "charge density (e/cm)"
plot 'pcdchar.dat' u 1:2:3 w errorbars 2, \
     'pcdchar.dat' u 1:2  notitle w histeps 1
pause -1

set title 'Path length distribution of transmitted particles'
set xlabel "path length (cm)"
set ylabel "PDF (1/cm)"
plot 'pctltr.dat' u 1:2:3 w errorbars 2,\
     'pctltr.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Path length distribution of backscattered particles'
set xlabel "path length (cm)"
set ylabel "PDF (1/cm)"
plot 'pctlbk.dat' u 1:2:3 w errorbars 2,\
     'pctlbk.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Path length distribution of absorbed particles'
set xlabel "path length (cm)"
set ylabel "PDF (1/cm)"
plot 'pctlab.dat' u 1:2:3 w errorbars 2,\
     'pctlab.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Angular distribution of emerging electrons'
set xlabel "angle (deg)"
set ylabel "PDF (1/sr)"
plot 'pcanel.dat' u 1:3:4 w errorbars 2,\
     'pcanel.dat' u 1:3 notitle w histeps 1
pause -1

set title 'Angular distribution of emerging photons'
set xlabel "angle (deg)"
set ylabel "PDF (1/sr)"
plot 'pcanga.dat' u 1:3:4 w errorbars 2,\
     'pcanga.dat' u 1:3 notitle w histeps 1
pause -1

set title 'Angular distribution of emerging positrons'
set xlabel "angle (deg)"
set ylabel "PDF (1/sr)"
plot 'pcanpo.dat' u 1:3:4 w errorbars 2,\
     'pcanpo.dat' u 1:3 notitle w histeps 1
pause -1

set title 'Energy distribution of transmitted electrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pceneltr.dat' u 1:2:3 w errorbars 2,\
     'pceneltr.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of backscattered electrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pcenelbk.dat' u 1:2:3 w errorbars 2,\
     'pcenelbk.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of transmitted photons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pcengatr.dat' u 1:2:3 w errorbars 2,\
     'pcengatr.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of backscattered photons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pcengabk.dat' u 1:2:3 w errorbars 2,\
     'pcengabk.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of transmitted positrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pcenpotr.dat' u 1:2:3 w errorbars 2,\
     'pcenpotr.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy distribution of backscattered positrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pcenpobk.dat' u 1:2:3 w errorbars 2,\
     'pcenpobk.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Deposited energy distribution in the first material'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
set logscale y
plot 'pcedepm1.dat' u 1:2:3 w errorbars 2,\
     'pcedepm1.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Deposited energy distribution in the first material'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
set nologscale y
plot 'pcedepm1.dat' u 1:2:3 w errorbars 2,\
     'pcedepm1.dat' u 1:2 notitle w histeps 1
pause -1

unset xzeroaxis

set style line 100 lt 5 lw 0.5
set pm3d solid hidden3d 100
set view 60, 150, 1,1
set ticslevel 0.07
set xlabel "z (cm)"
set ylabel "r (cm)"

set nologscale x
set nologscale y
set nologscale z
set xrange [*:*]
set yrange [*:*]
set zrange [*:*]

set mouse

set title '3D dose distribution (body 1)'
set zlabel "Dose (eV/g)"
splot 'pcdosch1.dat' u 1:2:3 w pm3d
pause -1 'Press enter to continue'

set title '3D deposited charge distribution (body 1)'
set zlabel "Charge (e/cm**3)"
splot 'pcdosch1.dat' u 1:2:5 w pm3d
pause -1 'Press enter to continue'

set title '3D dose distribution (body 2)'
set zlabel "Dose (eV/g)"
splot 'pcdosch2.dat' u 1:2:3 w pm3d
pause -1 'Press enter to continue'

set title '3D deposited charge distribution (body 2)'
set zlabel "Charge (e/cm**3)"
splot 'pcdosch2.dat' u 1:2:5 w pm3d
pause -1 'Press enter to continue'

set title '3D dose distribution (body 3)'
set zlabel "Dose (eV/g)"
splot 'pcdosch3.dat' u 1:2:3 w pm3d
pause -1 'Press enter to continue'

set title '3D deposited charge distribution (body 3)'
set zlabel "Charge (e/cm**3)"
splot 'pcdosch3.dat' u 1:2:5 w pm3d
pause -1 'Press enter to continue'

