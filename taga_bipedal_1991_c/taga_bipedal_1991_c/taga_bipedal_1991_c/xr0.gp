set datafile separator ","

set terminal postscript eps color enhanced
set out "xr0.eps"
set key

plot   "out.csv" using 1:2 title 'xr0' w l lw 4 lc rgb 'red',\
	   "out.csv" using 1:3 title 'xr'  w l lw 4,\
	   "out.csv" using 1:4 title 'yr0' w l lw 4 lc rgb 'red',\
	   "out.csv" using 1:5 title 'yr'  w l lw 4,\
	   "out.csv" using 1:6 title 'xl0' w l lw 4 lc rgb 'red',\
	   "out.csv" using 1:7 title 'xl'  w l lw 4,\
	   "out.csv" using 1:8 title 'yl0' w l lw 4 lc rgb 'red',\
	   "out.csv" using 1:9 title 'yl'  w l lw 4
