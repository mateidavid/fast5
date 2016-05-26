HDF_ROOT = /usr/local

.PHONY: all list clean

TARGETS = f5dump f5dump-full hdf5-mod f5-mod

all: ${TARGETS}

list:
	@echo "TARGETS=${TARGETS}"

clean:
	@rm -f ${TARGETS}

f5dump: f5dump.cpp fast5.hpp hdf5_tools.hpp
	g++ -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -pedantic -I ${HDF_ROOT}/include -o $@ $< -L ${HDF_ROOT}/lib -lhdf5 -Wl,--rpath=${HDF_ROOT}/lib

f5dump-full: f5dump-full.cpp fast5.hpp hdf5_tools.hpp
	g++ -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -pedantic -I ${HDF_ROOT}/include -o $@ $< -L ${HDF_ROOT}/lib -lhdf5 -Wl,--rpath=${HDF_ROOT}/lib

hdf5-mod: hdf5-mod.cpp fast5.hpp hdf5_tools.hpp
	g++ -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -pedantic -I ${HDF_ROOT}/include -o $@ $< -L ${HDF_ROOT}/lib -lhdf5 -Wl,--rpath=${HDF_ROOT}/lib

f5-mod: f5-mod.cpp fast5.hpp hdf5_tools.hpp
	g++ -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -pedantic -I ${HDF_ROOT}/include -o $@ $< -L ${HDF_ROOT}/lib -lhdf5 -Wl,--rpath=${HDF_ROOT}/lib
