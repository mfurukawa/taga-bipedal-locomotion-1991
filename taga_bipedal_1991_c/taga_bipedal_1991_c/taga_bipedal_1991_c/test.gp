set datafile separator ","

set terminal postscript eps
set out "img1.eps"
set key

plot   "out.csv" using 1:2 title 'x1' w l lw 4,\
	   "out.csv" using 1:3 title 'x2' w l lw 4,\
	   "out.csv" using 1:4 title 'x3' w l lw 4,\
	   "out.csv" using 1:5 title 'x4' w l lw 4,\
	   "out.csv" using 1:6 title 'x5' w l lw 8,\
	   "out.csv" using 1:7 title 'x6' w l lw 4,\
	   "out.csv" using 1:8 title 'x7' w l lw 4,\
	   "out.csv" using 1:9 title 'x8' w l lw 8,\
	   "out.csv" using 1:10 title 'x9' w l lw 4,\
	   "out.csv" using 1:11 title 'x10' w l lw 4,\
	   "out.csv" using 1:12 title 'u1' w l lw 4

