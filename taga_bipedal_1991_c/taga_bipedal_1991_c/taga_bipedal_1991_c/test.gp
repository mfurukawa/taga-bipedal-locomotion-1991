# Masahiro Furukawa
# m.furukawa@ist.osaka-u.ac.jp
#
# Oct 18, 2016
#
# taga1991
#
# reference : 
# http://www.ss.scphys.kyoto-u.ac.jp/person/yonezawa/contents/program/gnuplot/paper_adv1.html

set datafile separator ","

#set terminal postscript eps color enhanced
set terminal postscript eps
set out "img1.eps"
set key left

set size 1,1
set lmargin 10
set rmargin 2

set xlabel 'time'
set xrange [-0.1:10]
set yrange [-20:240]
set tics out

set nokey
set ylabel 'u(i)' 
set ytics 10
set ytics ('u_1' 0,\
		   'u_2' 20,\
		   'u_5' 40,\
		   'u_6' 60,\
		   'u_9' 80,\
		   'u_{10}' 100,\
		   'u_3' 120,\
		   'u_4' 140,\
		   'u_7' 160,\
		   'u_8' 180,\
		   'u_{11}' 200,\
		   'u_{12}' 220 )

plot   "out.csv" using 1:($24*1.1+0) title 'u1' w l lt -1 lw 4,\
	   "out.csv" using 1:($25*1.1+20) title 'u2' w l lt -1  lw 4,\
	   "out.csv" using 1:($28*1.1+40) title 'u5' w l lt -1  lw 4,\
	   "out.csv" using 1:($29*1.1+60) title 'u6' w l lt -1  lw 4,\
	   "out.csv" using 1:($32*1.1+80) title 'u9' w l lt -1  lw 4,\
	   "out.csv" using 1:($33*1.1+100) title 'u10' w l lt -1  lw 4,\
	   "out.csv" using 1:($26*1.1+120) title 'u3' w l lt -1  lw 4,\
	   "out.csv" using 1:($27*1.1+140) title 'u4' w l lt -1  lw 4,\
	   "out.csv" using 1:($30*1.1+160) title 'u7' w l lt -1  lw 4,\
	   "out.csv" using 1:($31*1.1+180) title 'u8' w l lt -1  lw 4,\
	   "out.csv" using 1:($34*1.1+200) title 'u11' w l lt -1  lw 4,\
	   "out.csv" using 1:($35*1.1+220) title 'u12' w l lt -1  lw 4

#set multiplot layout 3,1

#set ylabel 'x(i) position'

#plot   "out.csv" using 1:2 title 'x1' w l lw 4,\
#	   "out.csv" using 1:3 title 'x2' w l lw 4 dt(1,2) ,\
#	   "out.csv" using 1:4 title 'x3' w l lw 4,\
#	   "out.csv" using 1:5 title 'x4' w l lw 4 dt(1,2) ,\
#	   "out.csv" using 1:7 title 'x6' w l lw 4,\
#	   "out.csv" using 1:8 title 'x7' w l lw 4 dt(1,2) ,\
#	   "out.csv" using 1:9 title 'x8' w l lw 4,\
#	   "out.csv" using 1:11 title 'x10' w l lw 4 dt(1,2) ,\
#	   "out.csv" using 1:13 title 'x12' w l lw 4,\
#	   "out.csv" using 1:14 title 'x13' w l lw 4 dt(1,2)

#set tmargin 20
#set bmargin 00
#set ylabel 'x(i) as angle'

#plot   "out.csv" using 1:6 title 'x5' w l lw 4,\
#	   "out.csv" using 1:10 title 'x9' w l lw 4,\
#	   "out.csv" using 1:12 title 'x11' w l lw 4,\
#	   "out.csv" using 1:15 title 'x14' w l lw 4

#set tmargin 2

#unset multiplot