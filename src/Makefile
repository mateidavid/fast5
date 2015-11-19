HDF_ROOT = /usr/local

a: a.cpp fast5.hpp hdf5_tools.hpp
	g++ -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -pedantic -Wno-unused-parameter -I ${HDF_ROOT}/include -o $@ $^ -L ${HDF_ROOT}/lib -lhdf5 -Wl,--rpath=${HDF_ROOT}/lib
