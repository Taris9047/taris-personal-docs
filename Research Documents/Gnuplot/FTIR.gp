#!/usr/bin/env gnuplot --persist
set terminal postscript eps enhanced color "Helvetica" 24
#set terminal pdf enhanced fname 'Helvetica' fsize 12
#set terminal aqua
#set datafile separator ","
set style line 1 lw 4
set style line 2 lw 4

# Defining Margins
set lmargin 10
set bmargin 4
set rmargin 2
set tmargin 1

# Setting up Input and Output
IN_FILENAME = 'SiNxHadi.txt'
IN_FILENAME_BASELINE_FIT = 'Baseline_ary.txt'
IN_FILENAME_FTIR_SUB = 'FTIR_sub.txt'
OUT_FILENAME_BASELINE = 'FTIR_with_baseline.eps'
OUT_FILENAME_REFERENCE = 'FTIR_raw_data.eps'
OUT_FILENAME = 'FTIR_post_process.eps'

# Plotting the Raw Data
set output OUT_FILENAME_REFERENCE
set xlabel "Wavenumber (cm^{-1})"
set xtics 0,1000
set xrange [0:5500]
set ylabel "Absorbance (A.U.)"
set yrange [-0.8:0.6]
set format y ""
FTIR_TITLE = sprintf("FTIR Spectra of a-SiNx:H")
set title FTIR_TITLE
set key top right
show margin 
plot IN_FILENAME using ($1):($2) title "" w l 

# Plotting Baseline
set xlabel "Wavenumber (cm^{-1})"
set xtics 0,1000
set xrange [0:5500]
set ylabel "Absorbance (A.U.)"
set yrange [-0.8:0.6]
set format y ""
FTIR_TITLE = sprintf("FTIR Spectra of a-SiNx:H with Baseline")
set title FTIR_TITLE
set output OUT_FILENAME_BASELINE
set key top right
show margin
plot IN_FILENAME using ($1):($2) title "" w l, \
     IN_FILENAME_BASELINE_FIT using ($1):($2) title "Baseline" w l

# Plotting subtracted FTIR value
set xlabel "Wavenumber (cm^{-1})"
set xtics 0,1000
set xrange [0:5500]
set ylabel "Absorbance (A.U.)"
set yrange [-0.06:0.35]
set format y ""
FTIR_TITLE = sprintf("FTIR Spectra of a-SiNx:H (Baseline Process)")
set title FTIR_TITLE
set output OUT_FILENAME
set key top right
show margin
plot IN_FILENAME_FTIR_SUB using ($1):($2) title "" w l ls 1
