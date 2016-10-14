# Masahrio Furukawa
# m.furukawa@ist.osaka-u.ac.jp
#
# created Oct 14, 2016


CC = g++
LINK = -lm 
FLAG = -v -Wall

all: taga.o taga1991.o gauss_jordan.o
	$(CC) -o taga $^ $(LINK) 
gauss_jordan.o: gauss_jordan.cpp gauss_jordan.h
	$(CC) -c gauss_jordan.cpp -o $@ $(FLAG)
gauss_jordan_test.o: gauss_jordan_test.cpp gauss_jordan.cpp
	$(CC) -c $^ -o $@ $(FLAG)
taga1991.o: taga1991.cpp gauss_jordan.h
	$(CC) -o $@ -c taga1991.cpp $(FLAG)
taga.o: taga1991_test.cpp
	$(CC) -o $@ -c $^ $(FLAG)

clean:
	rm *.o
	rm taga
