#  -------------------
#  Usage: 
#         gmake clean
#         gmake 
#  -------------------
#  GNU_Makefile to 
#  - COMPILE and LINK  test_mocsy.f90 and mocsy.f90
#  - generate library libmocsy.a
#  - generate python interface module _mocsy.so
#
#  James Orr, LSCE/IPSL, CEA-CNRS-UVSQ, 91191 Gif-sur-Yvette, France
#  15 January 2014

#=======================================================================
#                   define the compiler names
#=======================================================================
CC       = gcc
# To use another Fortran compiler, replace "f95" in FC and F90 lines with your compiler command
# For example, comment out 2 lines below with "f95" & uncomment the following 2 lines for "gfortran"
#FC = fort77
#FC = xlf
#FC = f95
#F90 = f95
F90 = gfortran
#FC = ifort
#F90 = ifort
PYTHON   = python2

#DEBUGFLAGS = -g
#LDFLAGS = -L/usr/local/lib -lnetcdf -lnetcdff
LDFLAGS = -L./ -lmocsy
INCLUDEFLAGS = -I/usr/local/include .


#=======================================================================
#                   define desired precision
#=======================================================================

# set to  2 if you wish results in DOUBLE precision
# set to 1 or 0 if SINGLE
PRECISION = 2

#=======================================================================
#                     additional flags
#=======================================================================

ifeq ($(F90),gfortran)
	FPP      = gfortran -E
	FPP_F90FLAGS = -x f95-cpp-input -fPIC -DUSE_PRECISION=$(PRECISION)
	F90FLAGS = -fPIC -x f95-cpp-input -DUSE_PRECISION=$(PRECISION)
    #FCOMP    = gfortran
    FCOMP    = gnu95
    LIBS     =
endif

ifeq ($(F90),ifort)

	FPP      = gfortran -E # gfortran f90wrap temp files only. not compilation
	FPP_F90FLAGS = -x f95-cpp-input -fPIC -DUSE_PRECISION=$(PRECISION)
	F90FLAGS = -fpscomp logicals -fPIC -DUSE_PRECISION=$(PRECISION)  # use 1 and 0 for True and False
    FCOMP    = intelem # for f2py
    LIBS =
endif

CFLAGS = -fPIC #     ==> universal for ifort, gfortran, pgi

#=======================================================================
#=======================================================================

UNAME = $(shell uname)

ifeq (${UNAME}, Darwin)
  LIBTOOL = libtool -static -o
else
  LIBTOOL = ar src
endif

# ======================================================================
# PROJECT CONFIG, do not put spaced behind the variables
# ======================================================================
# mapping between Fortran and C types
KIND_MAP = kind_map

#=======================================================================
#       List all source files required for the project
#=======================================================================

LIBSRC_SOURCES = \
          singledouble \
	  DNAD \
	  sw_adtg \
          sw_ptmp \
          sw_temp \
          tpot \
          tis \
          p80 \
          phsolvers \
          rho \
          rhoinsitu \
          depth2press \
          constants \
          varsolver \
          vars \
          p2fCO2 \
          f2pCO2 \
          gasx \
	  derivnum

# file names
LIBSRC_FILES = $(addsuffix .f90,${LIBSRC_SOURCES})

# object files
LIBSRC_OBJECTS = $(addsuffix .o,${LIBSRC_SOURCES})

# only used when cleaning up
LIBSRC_FPP_FILES = $(addsuffix .fpp,${LIBSRC_SOURCES})

#=======================================================================
#       List all source files that require a Python interface
#=======================================================================

# names (without suffix), f90 sources
LIBSRC_WRAP_SOURCES = ${LIBSRC_SOURCES}

# file names
LIBSRC_WRAP_FILES = $(addsuffix .f90,${LIBSRC_WRAP_SOURCES})

# object files
LIBSRC_WRAP_OBJECTS = $(addsuffix .o,${LIBSRC_WRAP_SOURCES})

# fpp files
LIBSRC_WRAP_FPP_FILES = $(addsuffix .fpp,${LIBSRC_WRAP_SOURCES})


#=======================================================================
#
#=======================================================================

# List of executables to be built within the package
PROGRAMS = libmocsy.a test_mocsy test_vars _mocsy.so

# "make" builds all
all: $(PROGRAMS)

#---------------------------------------------------------------------------
EXEC = test_mocsy

library = libmocsy.a

#---------------------------------------------------------------------------
# Build the mocsy library containing the object files (not used, illustration only)
$(library):  $(LIBSRC_OBJECTS)
	ar cr $(library) $(LIBSRC_OBJECTS)

# Build the Fortran program executable that tests the mocsy library (test_mocsy)
$(EXEC): $(LIBSRC_OBJECTS) test_mocsy.o $(library) 
	${F90} ${F90FLAGS} -o $@ $@.f90 $(LDFLAGS)

test_vars:  $(LIBSRC_OBJECTS) test_vars.o $(library) 
	${F90} ${F90FLAGS} -o $@ $@.f90 $(LDFLAGS)

test_derivauto:  $(LIBSRC_OBJECTS) test_derivauto.o $(library) 
	${F90} ${F90FLAGS} -o $@ $@.f90 $(LDFLAGS)

test_derivnum:  $(LIBSRC_OBJECTS) test_derivnum.o $(library) 
	${F90} ${F90FLAGS} -o $@ $@.f90 $(LDFLAGS)

#=======================================================================
#                 Relevant suffixes
#=======================================================================

.SUFFIXES: .f90 .fpp

#---------------------------------------------------------------------------

# Utility targets
.PHONY: clean veryclean


clean:
	-rm *.o *.so *.a *.mod *.fpp f90wrap*.f90 mocsy.py *.pyc
	-rm -rf mocsy_pkg/ src.*


.f90.o:
	${F90} ${F90FLAGS} -c $< -o $@


.c.o:
	${CC} ${CFLAGS} -c $< -o $@


.f90.fpp:
	${FPP} ${FPP_F90FLAGS} $<  -o $@


libsrc.a: ${LIBSRC_OBJECTS}
	${LIBTOOL} $@ $?


_mocsy.so: libsrc.a ${LIBSRC_FPP_FILES}
	f90wrap -m mocsy ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -v
	f2py-f90wrap --fcompiler=$(FCOMP) --build-dir . -c -m _mocsy -L. -lsrc f90wrap*.f90


_mocsy_pkg.so: libsrc.a ${LIBSRC_FPP_FILES}
	f90wrap -m mocsy_pkg ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -v -P
	f2py-f90wrap --fcompiler=$(FCOMP) --build-dir . -c -m _mocsy_pkg -L. -lsrc f90wrap*.f90

veryclean: clean
	rm -f *~ $(PROGRAMS) 
