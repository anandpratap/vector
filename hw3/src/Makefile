serial:
	icpc -O3 -vec-report3 -o serial.o main_serial.cpp
parallel:
	mpicxx -O0 -Wall -vec-report3 -o parallel.o main.cpp
debug:
	mpicxx -O0 -vec-report3 -g -debug -o parallel.o main.cpp 
clean:
	rm *.o*
