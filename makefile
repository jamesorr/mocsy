#  -------------------
#  Usage: 
#         gmake clean
#         gmake 
#  -------------------
#  GNU_Makefile to COMPILE and LINK  test_mocsy.f90 and mocsy.f90
#  James Orr, LSCE/IPSL, CEA-CNRS-UVSQ, 91191 Gif-sur-Yvette, France
#  15 January 2014

# To use another Fortran compiler, replace "f95" in FC and F90 lines with your compiler command
# For example, comment out 2 lines below with "f95" & uncomment the following 2 lines for "gfortran"
#FC = fort77
#FC = xlf
#FC = f95
#F90 = f95
FC = gfortran
F90 = gfortran
#FC = ifort
#F90 = ifort

#DEBUGFLAGS = -g
#LDFLAGS = -L/usr/local/lib -lnetcdf -lnetcdff
LDFLAGS = -L./ -lmocsy
INCLUDEFLAGS = -I/usr/local/include .

# List of executables to be built within the package
PROGRAMS = libmocsy.a test_mocsy test_mocsy2 mocsy.so

# "make" builds all
all: $(PROGRAMS)

#---------------------------------------------------------------------------

SOURCES = singledouble.f90 \
          sw_adtg.f90 \
          sw_ptmp.f90 \
          sw_temp.f90 \
          tpot.f90 \
          tis.f90 \
          p80.f90 \
          phsolvers.f90 \
          rho.f90 \
          rhoinsitu.f90 \
          depth2press.f90 \
          constants.f90 \
          varsolver.f90 \
          vars_rho1028.f90 \
          vars_rho1000.f90 \
          vars.f90 \
          vars2.f90 \
          p2fCO2.f90 \
          f2pCO2.f90 \
          gasx.f90 

OBJS =  singledouble.o \
        sw_adtg.o \
        sw_ptmp.o \
        sw_temp.o \
        tpot.o \
        tis.o \
        p80.o \
        phsolvers.o \
        rho.o \
        rhoinsitu.o \
        depth2press.o \
        constants.o \
        varsolver.o \
        vars_rho1028.o \
        vars_rho1000.o \
        vars.o \
        vars2.o \
        p2fCO2.o \
        f2pCO2.o \
        gasx.o

#TOBJS = $(OBJS) \
#        test_mocsy.o

#TOBJS2 = $(OBJS) \
#        test_mocsy2.o

EXEC = test_mocsy

EXEC2 = test_mocsy2

library = libmocsy.a
#---------------------------------------------------------------------------

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
#%: %.o
#	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General Pattern rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

#---------------------------------------------------------------------------
# Build the mocsy library containing the object files (not used, illustration only)
$(library):  $(OBJS)
	ar cr $(library) $(OBJS)

# Build the Fortran program executable that tests the mocsy library (test_mocsy)
#$(EXEC): $(TOBJS) $(library)
#	$(FC) $(FCFLAGS) -o $(EXEC) $(TOBJS) 
$(EXEC): $(OBJS) $(library) 
	$(FC) $(FCFLAGS) -o $@ $@.f90 $(LDFLAGS)

# Build the Fortran program executable that tests the mocsy library (test_mocsy)
#$(EXEC2): $(TOBJS2) $(library)
#	$(FC) $(FCFLAGS) -o $(EXEC2) $(TOBJS2) 
#$(EXEC2): $(OBJS) $(library) $(EXEC2).o
#	$(FC) $(FCFLAGS) -o $(EXEC2) $(EXEC2).o $(LDFLAGS)
$(EXEC2): $(OBJS) $(library) 
	$(FC) $(FCFLAGS) -o $@ $@.f90 $(LDFLAGS)

# Build the shared object file for python
mocsy.so: $(OBJS)
	f2py -c $(SOURCES) -m mocsy --fcompiler=gnu95 --f90flags=-O3
#---------------------------------------------------------------------------

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ 

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.so

veryclean: clean
	rm -f *~ $(PROGRAMS) 
