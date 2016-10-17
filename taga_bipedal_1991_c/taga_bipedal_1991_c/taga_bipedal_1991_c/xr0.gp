set datafile separator ","

set terminal postscript eps color enhanced
set out "xr0.eps"
set key

plot   "out.csv" using 1:2 title 'xr0' w l lw 4 lc rgb 'red',\
	   "out.csv" using 1:3 title 'xr' w l lw 4,\
	   "out.csv" using 1:4 title 'yr0' w l lw 4,\
	   "out.csv" using 1:5 title 'yr' w l lw 4
