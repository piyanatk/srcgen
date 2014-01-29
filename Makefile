############################################
# Makefile for srcgen
# Created: 2014-01-23 by Piyanat Kittiwisit
############################################


# ==========================================
# Required libraries
# ------------------------------------------
# You might have to change the following lines.
# srcgen requires HEALPix <= version 2.15a. Newer version will not work due
# change in some libraries.You must build the c++ library before compiling simsky
HPX_CXX_INC = ./Healpix_2.15a/src/cxx/generic_gcc/include
HPX_CXX_LIB = ./Healpix_2.15a/src/cxx/generic_gcc/lib
# ==========================================cd

CXX = g++
CXX_INCS = -I$(HPX_CXX_INC) -L$(HPX_CXX_LIB)
CXXFLAGS = -g -Wall $(CXX_INCS)

all: srcgen

srcgen:
	$(CXX) $(CXXFLAGS) -o srcgen src/srcgen.cc -lgsl -lgslcblas -lm -lhealpix_cxx -lcxxsupport -lfftpack -lcfitsio

clean:
	rm -rf src/*.o srcgen
