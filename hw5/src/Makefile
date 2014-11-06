all:
	icpc -O0 -o serial serial.cpp
	icpc -O0 -openmp -O3 -o parallel parallel.cpp
debug:
	icpc -g -debug -o serial serial.cpp
