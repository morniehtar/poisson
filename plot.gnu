#! /usr/bin/gnuplot

set encoding utf8
set tmargin 1
set tics out
set grid
set term qt title "Potential distribution" font "arial,12" enhanced #persist

set style line 1 lt rgb "#ff0000" lw 2 #Scarlet, thick
set style line 2 lt rgb "#ff8000" lw 2 #Ginger, thick
set style line 3 lt rgb "#ffbf00" lw 2 #Orange, thick
set style line 4 lt rgb "#ffff00" lw 2 #Yellow, thick
set style line 5 lt rgb "#80ff00" lw 2 #Lime, thick
set style line 6 lc rgb "#00ffff" lw 1.2 dt 2 #Turquoise, thin, dotted

plot "./output.dat" using 1:2 smooth csplines linestyle 1 title "0", \
"./output.dat" using 1:3 smooth csplines linestyle 2 title "{/symbol p}/6", \
"./output.dat" using 1:4 smooth csplines linestyle 3 title "{/symbol p}/4", \
"./output.dat" using 1:5 smooth csplines linestyle 4 title "{/symbol p}/3", \
"./output.dat" using 1:6 smooth csplines linestyle 5 title "{/symbol p}/2", \
"./output.dat" using 1:7 smooth csplines linestyle 6 title "{/symbol p}"

pause mouse close
