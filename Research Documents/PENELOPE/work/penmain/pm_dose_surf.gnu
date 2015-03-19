# Gnuplot MS-Windows 32 bit version 4.0
# Plots results from 'penmain.f'

# unset mouse

set title 'Depth-dose distribution'
set nologscale y
set xrange [*:*]
set yrange [0.:*]
set ylabel 'Depth-dose (eV/(g/cm**2)) '
set xlabel 'z (cm)'
plot 'pm_depth_dose.dat' u 1:2:3 w errorbars 2,\
     'pm_depth_dose.dat' u 1:2 notitle w histeps 1
pause -1

set zero 1.0e-60
set style line 100 lt 5 lw 0.5
set pm3d solid hidden3d 100
set view 50, 60, 1,1
set ticslevel 0.05
set xlabel "x (cm)"
set ylabel "y (cm)"
set zlabel "Dose (eV/g)"

set nologscale x
set nologscale y
set nologscale z
set xrange [*:*]
set yrange [*:*]
set zrange [0.:*]

set title '2D dose distribution (ev/g). Plane I3=1' ; splot 'pm_2d_dose_1.dat'  u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=2' ; splot 'pm_2d_dose_2.dat'  u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=3' ; splot 'pm_2d_dose_3.dat'  u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=4' ; splot 'pm_2d_dose_4.dat'  u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=5' ; splot 'pm_2d_dose_5.dat'  u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=6' ; splot 'pm_2d_dose_6.dat'  u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=7' ; splot 'pm_2d_dose_7.dat'  u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=8' ; splot 'pm_2d_dose_8.dat'  u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=9' ; splot 'pm_2d_dose_9.dat'  u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=10'; splot 'pm_2d_dose_10.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=11'; splot 'pm_2d_dose_11.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=12'; splot 'pm_2d_dose_12.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=13'; splot 'pm_2d_dose_13.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=14'; splot 'pm_2d_dose_14.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=15'; splot 'pm_2d_dose_15.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=16'; splot 'pm_2d_dose_16.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=17'; splot 'pm_2d_dose_17.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=18'; splot 'pm_2d_dose_18.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=19'; splot 'pm_2d_dose_19.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=20'; splot 'pm_2d_dose_20.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=21'; splot 'pm_2d_dose_21.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=22'; splot 'pm_2d_dose_22.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=23'; splot 'pm_2d_dose_23.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=24'; splot 'pm_2d_dose_24.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=25'; splot 'pm_2d_dose_25.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=26'; splot 'pm_2d_dose_26.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=27'; splot 'pm_2d_dose_27.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=28'; splot 'pm_2d_dose_28.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=29'; splot 'pm_2d_dose_29.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=30'; splot 'pm_2d_dose_30.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=31'; splot 'pm_2d_dose_31.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=32'; splot 'pm_2d_dose_32.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=33'; splot 'pm_2d_dose_33.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=34'; splot 'pm_2d_dose_34.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=35'; splot 'pm_2d_dose_35.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=36'; splot 'pm_2d_dose_36.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=37'; splot 'pm_2d_dose_37.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=38'; splot 'pm_2d_dose_38.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=39'; splot 'pm_2d_dose_39.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=40'; splot 'pm_2d_dose_40.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=41'; splot 'pm_2d_dose_41.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=42'; splot 'pm_2d_dose_42.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=43'; splot 'pm_2d_dose_43.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=44'; splot 'pm_2d_dose_44.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=45'; splot 'pm_2d_dose_45.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=46'; splot 'pm_2d_dose_46.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=47'; splot 'pm_2d_dose_47.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=48'; splot 'pm_2d_dose_48.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=49'; splot 'pm_2d_dose_49.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=50'; splot 'pm_2d_dose_50.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=51'; splot 'pm_2d_dose_51.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=52'; splot 'pm_2d_dose_52.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=53'; splot 'pm_2d_dose_53.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=54'; splot 'pm_2d_dose_54.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=55'; splot 'pm_2d_dose_55.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=56'; splot 'pm_2d_dose_56.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=57'; splot 'pm_2d_dose_57.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=58'; splot 'pm_2d_dose_58.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=59'; splot 'pm_2d_dose_59.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=60'; splot 'pm_2d_dose_60.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=61'; splot 'pm_2d_dose_61.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=62'; splot 'pm_2d_dose_62.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=63'; splot 'pm_2d_dose_63.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=64'; splot 'pm_2d_dose_64.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=65'; splot 'pm_2d_dose_65.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=66'; splot 'pm_2d_dose_66.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=67'; splot 'pm_2d_dose_67.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=68'; splot 'pm_2d_dose_68.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=69'; splot 'pm_2d_dose_69.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=70'; splot 'pm_2d_dose_70.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=71'; splot 'pm_2d_dose_71.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=72'; splot 'pm_2d_dose_72.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=73'; splot 'pm_2d_dose_73.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=74'; splot 'pm_2d_dose_74.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=75'; splot 'pm_2d_dose_75.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=76'; splot 'pm_2d_dose_76.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=77'; splot 'pm_2d_dose_77.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=78'; splot 'pm_2d_dose_78.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=79'; splot 'pm_2d_dose_79.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=80'; splot 'pm_2d_dose_80.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=81'; splot 'pm_2d_dose_81.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=82'; splot 'pm_2d_dose_82.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=83'; splot 'pm_2d_dose_83.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=84'; splot 'pm_2d_dose_84.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=85'; splot 'pm_2d_dose_85.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=86'; splot 'pm_2d_dose_86.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=87'; splot 'pm_2d_dose_87.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=88'; splot 'pm_2d_dose_88.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=89'; splot 'pm_2d_dose_89.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=90'; splot 'pm_2d_dose_90.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=91'; splot 'pm_2d_dose_91.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=92'; splot 'pm_2d_dose_92.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=93'; splot 'pm_2d_dose_93.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=94'; splot 'pm_2d_dose_94.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=95'; splot 'pm_2d_dose_95.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=96'; splot 'pm_2d_dose_96.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=97'; splot 'pm_2d_dose_97.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=98'; splot 'pm_2d_dose_98.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=99'; splot 'pm_2d_dose_99.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=100'; splot 'pm_2d_dose_100.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=101'; splot 'pm_2d_dose_101.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=102'; splot 'pm_2d_dose_102.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=103'; splot 'pm_2d_dose_103.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=104'; splot 'pm_2d_dose_104.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=105'; splot 'pm_2d_dose_105.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=106'; splot 'pm_2d_dose_106.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=107'; splot 'pm_2d_dose_107.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=108'; splot 'pm_2d_dose_108.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=109'; splot 'pm_2d_dose_109.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=110'; splot 'pm_2d_dose_110.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=111'; splot 'pm_2d_dose_111.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=112'; splot 'pm_2d_dose_112.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=113'; splot 'pm_2d_dose_113.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=114'; splot 'pm_2d_dose_114.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=115'; splot 'pm_2d_dose_115.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=116'; splot 'pm_2d_dose_116.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=117'; splot 'pm_2d_dose_117.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=118'; splot 'pm_2d_dose_118.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=119'; splot 'pm_2d_dose_119.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

set title '2D dose distribution (eV/g). Plane I3=120'; splot 'pm_2d_dose_120.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=121'; splot 'pm_2d_dose_121.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=122'; splot 'pm_2d_dose_122.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=123'; splot 'pm_2d_dose_123.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=124'; splot 'pm_2d_dose_124.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=125'; splot 'pm_2d_dose_125.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=126'; splot 'pm_2d_dose_126.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=127'; splot 'pm_2d_dose_127.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=128'; splot 'pm_2d_dose_128.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'
set title '2D dose distribution (eV/g). Plane I3=129'; splot 'pm_2d_dose_129.dat' u 4:5:7 notitle w pm3d; pause -1 'Press enter to continue'

