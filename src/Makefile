.SUFFIXES:
MAKEFLAGS += -r
SHELL := /bin/bash
.DELETE_ON_ERROR:
.PHONY: all help list clean

HDF5_DIR = /usr/local
HDF5_INCLUDE_DIR = ${HDF5_DIR}/include
HDF5_LIB_DIR = ${HDF5_DIR}/lib
HDF5_LIB = hdf5

TARGETS = f5dump f5dump-full hdf5-mod f5-mod

all: ${TARGETS}

print-%:
	@echo '$*=$($*)'

help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

list:
	@echo "TARGETS=${TARGETS}"

clean:
	@rm -f ${TARGETS}

check_hdf5:
	@[ -f "${HDF5_INCLUDE_DIR}/H5pubconf.h" ] || { echo "HDF5 headers not found" >&2; exit 1; }
	@[ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.so" ] || [ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.a" ] || { echo "HDF5 library not found" >&2; exit 1; }

%: %.cpp fast5.hpp hdf5_tools.hpp check_hdf5
	g++ -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -pedantic -I${HDF5_INCLUDE_DIR} -o $@ $< -L${HDF5_LIB_DIR} -Wl,--rpath=${HDF5_LIB_DIR} -l${HDF5_LIB} -lpthread -lz -ldl
