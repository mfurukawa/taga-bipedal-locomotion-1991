# Masahrio Furukawa
# m.furukawa@ist.osaka-u.ac.jp
#
# created Oct 14, 2016


CC = g++
LINK = -lm 
FLAG = -v

all: taga.o taga1991.o gauss_jordan.o
	$(CC) -o taga $^ $(LINK) 
gauss_jordan.o: gauss_jordan.cpp
	$(CC) -c $^ -o $@ $(FLAG)
gauss_jordan_test.o: gauss_jordan_test.cpp
	$(CC) -c $^ -o $@ $(FLAG)
taga1991.o: taga1991.cpp
	$(CC) -o $@ -c $^ $(FLAG)
taga.o: taga1991_test.cpp
	$(CC) -o $@ -c $^ $(FLAG)

clean:
	rm *.o
	rm taga
