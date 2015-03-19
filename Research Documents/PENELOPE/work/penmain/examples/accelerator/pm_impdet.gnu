# Gnuplot MS-Windows 32 bit version 4.0
# Plots results from 'penmain.f'

# unset mouse

set title 'Energy spectrum from impulse detector #1'
set nologscale y
set xrange [*:*]
set yrange [*:*]
set ylabel 'P(E)'
set xlabel 'E'
plot 'pm_spc_impdet_1.dat' u 1:2:3 w errorbars 2,\
     'pm_spc_impdet_1.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy spectrum from impulse detector #1'
set logscale y
set xrange [*:*]
set yrange [1.0e-12:*]
set ylabel 'P(E)'
set xlabel 'E'
plot 'pm_spc_impdet_1.dat' u 1:2:3 w errorbars 2,\
     'pm_spc_impdet_1.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy spectrum from impulse detector #2'
set nologscale y
set xrange [*:*]
set yrange [*:*]
set ylabel 'P(E)'
set xlabel 'E'
plot 'pm_spc_impdet_2.dat' u 1:2:3 w errorbars 2,\
     'pm_spc_impdet_2.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy spectrum from impulse detector #2'
set logscale y
set xrange [*:*]
set yrange [1.0e-12:*]
set ylabel 'P(E)'
set xlabel 'E'
plot 'pm_spc_impdet_2.dat' u 1:2:3 w errorbars 2,\
     'pm_spc_impdet_2.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy spectrum from impulse detector #3'
set nologscale y
set xrange [*:*]
set yrange [*:*]
set ylabel 'P(E)'
set xlabel 'E'
plot 'pm_spc_impdet_3.dat' u 1:2:3 w errorbars 2,\
     'pm_spc_impdet_3.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy spectrum from impulse detector #3'
set logscale y
set xrange [*:*]
set yrange [1.0e-12:*]
set ylabel 'P(E)'
set xlabel 'E'
plot 'pm_spc_impdet_3.dat' u 1:2:3 w errorbars 2,\
     'pm_spc_impdet_3.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy spectrum from impulse detector #4'
set nologscale y
set xrange [*:*]
set yrange [*:*]
set ylabel 'P(E)'
set xlabel 'E'
plot 'pm_spc_impdet_4.dat' u 1:2:3 w errorbars 2,\
     'pm_spc_impdet_4.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy spectrum from impulse detector #4'
set logscale y
set xrange [*:*]
set yrange [1.0e-12:*]
set ylabel 'P(E)'
set xlabel 'E'
plot 'pm_spc_impdet_4.dat' u 1:2:3 w errorbars 2,\
     'pm_spc_impdet_4.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy spectrum from impulse detector #5'
set nologscale y
set xrange [*:*]
set yrange [*:*]
set ylabel 'P(E)'
set xlabel 'E'
plot 'pm_spc_impdet_5.dat' u 1:2:3 w errorbars 2,\
     'pm_spc_impdet_5.dat' u 1:2 notitle w histeps 1
pause -1

set title 'Energy spectrum from impulse detector #5'
set logscale y
set xrange [*:*]
set yrange [1.0e-12:*]
set ylabel 'P(E)'
set xlabel 'E'
plot 'pm_spc_impdet_5.dat' u 1:2:3 w errorbars 2,\
     'pm_spc_impdet_5.dat' u 1:2 notitle w histeps 1
pause -1

