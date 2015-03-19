# Gnuplot MS-Windows 32 bit version 4.0
# Plots results from 'penmain.f'

# unset mouse
set zero 1.0e-60

set xzeroaxis
set yrange [:]

set title 'Energy distribution of transmitted electrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pm_energy_el_trans.dat' u 1:2:3 w errorbars 2,\
     'pm_energy_el_trans.dat' u 1:2 notitle w histeps 1
pause -1 'Press enter to continue'

set title 'Energy distribution of backscattered electrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pm_energy_el_back.dat' u 1:2:3 w errorbars 2,\
     'pm_energy_el_back.dat' u 1:2 notitle w histeps 1
pause -1 'Press enter to continue'

set title 'Energy distribution of transmitted photons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pm_energy_ph_trans.dat' u 1:2:3 w errorbars 2,\
     'pm_energy_ph_trans.dat' u 1:2 notitle w histeps 1
pause -1 'Press enter to continue'

set title 'Energy distribution of backscattered photons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pm_energy_ph_back.dat' u 1:2:3 w errorbars 2,\
     'pm_energy_ph_back.dat' u 1:2 notitle w histeps 1
pause -1 'Press enter to continue'

set title 'Energy distribution of transmitted positrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pm_energy_po_trans.dat' u 1:2:3 w errorbars 2,\
     'pm_energy_po_trans.dat' u 1:2 notitle w histeps 1
pause -1 'Press enter to continue'

set title 'Energy distribution of backscattered positrons'
set xlabel "energy (eV)"
set ylabel "PDF (1/eV)"
plot 'pm_energy_po_back.dat' u 1:2:3 w errorbars 2,\
     'pm_energy_po_back.dat' u 1:2 notitle w histeps 1
pause -1 'Press enter to continue'

unset xzeroaxis

set style line 100 lt 5 lw 0.5
set pm3d solid hidden3d 100
set view 50, 155, 1,1
set ticslevel 0.07
set xlabel "theta (deg)"
set ylabel "phi (deg)"
set zlabel "PDF (1/sr)"

set nologscale x
set nologscale y
set nologscale z
set xrange [1:180]
set yrange [0:360]
set zrange [*:*]

set mouse

set title 'Angular distribution of emerging electrons'
splot 'pm_angle_el.dat' u 1:2:3 notitle w pm3d,'pm_angle_el.dat' u 1:2:3 w l 1
pause -1 'Press enter to continue'

set title 'Angular distribution of emerging photons'
splot 'pm_angle_ph.dat' u 1:2:3 notitle w pm3d,'pm_angle_ph.dat' u 1:2:3 w l 1
pause -1 'Press enter to continue'

set title 'Angular distribution of emerging positrons'
splot 'pm_angle_po.dat' u 1:2:3 notitle w pm3d,'pm_angle_po.dat' u 1:2:3 w l 1
pause -1 'Press enter to continue'

