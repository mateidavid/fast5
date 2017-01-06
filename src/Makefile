.SUFFIXES:
MAKEFLAGS += -r
SHELL := /bin/bash
.DELETE_ON_ERROR:
.PHONY: all help list clean check_hdf5

HDF5_DIR ?= /usr/local
HDF5_INCLUDE_DIR ?= ${HDF5_DIR}/include
HDF5_LIB_DIR ?= ${HDF5_DIR}/lib
HDF5_LIB ?= hdf5
TCLAP_DIR ?= tclap
HPPTOOLS_DIR ?= hpptools

TARGETS = f5ls f5ls-full hdf5-mod f5-mod

all: ${TARGETS}

print-%:
	@echo '$*=$($*)'

help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

list:
	@echo "TARGETS=${TARGETS}"

clean:
	rm -rf ${TARGETS}

check_hdf5:
	@[ -f "${HDF5_INCLUDE_DIR}/H5pubconf.h" ] || { echo "HDF5 headers not found" >&2; exit 1; }
	@[ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.so" ] || [ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.a" ] || { echo "HDF5 library not found" >&2; exit 1; }

check_tclap:
	@[ -f "${TCLAP_DIR}/include/tclap/CmdLine.h" ] || { echo "TCLAP not found; get it from https://github.com/mateidavid/tclap.git" >&2; exit 1; }

check_hpptools:
	@[ -f "${HPPTOOLS_DIR}/include/alg.hpp" ] || { echo "HPPTOOLS not found; get it from https://github.com/mateidavid/hpptools.git" >&2; exit 1; }

%: %.cpp fast5.hpp hdf5_tools.hpp | check_hdf5
	${CXX} -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -Wpedantic -isystem ${HDF5_INCLUDE_DIR} -o $@ $< -L${HDF5_LIB_DIR} -Wl,--rpath=${HDF5_LIB_DIR} -l${HDF5_LIB} -lpthread -lz -ldl

f5dump: f5dump.cpp fast5.hpp hdf5_tools.hpp | check_hdf5 check_tclap check_hpptools
	${CXX} -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -Wpedantic -isystem ${HDF5_INCLUDE_DIR} -isystem ${TCLAP_DIR}/include -I ${HPPTOOLS_DIR}/include -o $@ $< -L${HDF5_LIB_DIR} -Wl,--rpath=${HDF5_LIB_DIR} -l${HDF5_LIB} -lpthread -lz -ldl

f5pack: f5pack.cpp fast5.hpp hdf5_tools.hpp | check_hdf5 check_tclap check_hpptools
	${CXX} -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -Wpedantic -isystem ${HDF5_INCLUDE_DIR} -isystem ${TCLAP_DIR}/include -I ${HPPTOOLS_DIR}/include -o $@ $< -L${HDF5_LIB_DIR} -Wl,--rpath=${HDF5_LIB_DIR} -l${HDF5_LIB} -lpthread -lz -ldl
