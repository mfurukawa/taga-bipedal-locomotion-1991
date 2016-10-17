set datafile separator ","

set terminal postscript eps color enhanced
set out "sp.eps"
set key

plot  "out.csv" using 1:2 title 'sp1' w l lw 4 lc rgb 'red'
