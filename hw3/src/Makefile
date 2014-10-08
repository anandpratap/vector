serial:
	icpc -O3 -vec-report3 -o serial.o main_serial.cpp
parallel:
	mpiicpc -O0 -Wall -vec-report3 -o parallel.o main.cpp
debug:
	mpiicpc -O0 -vec-report3 -g -debug -o parallel.o main.cpp 
clean:
	rm *.o*
