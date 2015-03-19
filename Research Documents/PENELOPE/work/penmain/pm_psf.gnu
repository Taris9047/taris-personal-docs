# Gnuplot MS-Windows 32 bit version 4.0
# Plots results from 'penmain.f'

unset mouse
set zero 1.0e-60
set size ratio -1

set title 'Phase-space file, detector #1, x-y distribution'
set nologscale x
set nologscale y
set xrange [*:*]
set yrange [*:*]
set ylabel 'x (cm)'
set xlabel 'y (cm)'
plot 'pm_psf_impdet_1.dat' u 3:4 w d 1
pause -1

set title 'Phase-space file, detector #1, x-y-E distribution'
set ticslevel 0.1
set nologscale x
set nologscale y
set xrange [*:*]
set yrange [*:*]
set ylabel 'x (cm)'
set xlabel 'y (cm)'
set zlabel 'E (eV)'
splot 'pm_psf_impdet_1.dat' u 3:4:2 w d 1
pause -1

set title 'Phase-space file, detector #2, x-y distribution'
set nologscale x
set nologscale y
set xrange [*:*]
set yrange [*:*]
set ylabel 'x (cm)'
set xlabel 'y (cm)'
plot 'pm_psf_impdet_2.dat' u 3:4 w d 1
pause -1

set title 'Phase-space file, detector #2, x-y-E distribution'
set ticslevel 0.1
set nologscale x
set nologscale y
set xrange [*:*]
set yrange [*:*]
set ylabel 'x (cm)'
set xlabel 'y (cm)'
set zlabel 'E (eV)'
splot 'pm_psf_impdet_2.dat' u 3:4:2 w d 1
pause -1

