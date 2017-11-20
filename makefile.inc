#$Id: makefile_itp,v 1.1.1.1 2016/05/30 06:22:14 zjcao Exp $
filein  =  

LDLIBS  = -lifcore

CXXAPPFLAGS = -O2 -Wno-deprecated -Dfortran3 -Dnewc
f90appflags = -O2 -fpp
f90 = mpif90
f77 = mpif77
CXX = mpicxx
CC  = mpicc
CLINKER = mpicxx 
