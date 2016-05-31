HDF5_ROOT = /usr/local

.PHONY: all list clean

TARGETS = f5dump f5dump-full hdf5-mod f5-mod

all: ${TARGETS}

list:
	@echo "TARGETS=${TARGETS}"

clean:
	@rm -f ${TARGETS}

%: %.cpp fast5.hpp hdf5_tools.hpp
	g++ -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -pedantic -I ${HDF5_ROOT}/include -o $@ $< -L ${HDF5_ROOT}/lib -Wl,--rpath=${HDF5_ROOT}/lib -lhdf5 -lpthread -lz -ldl
