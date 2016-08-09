#  -------------------
#  Usage: 
#         make clean
#         make 
#  -------------------
#  GNU_Makefile to COMPILE and LINK  test_mocsy.f90 and mocsy.f90
#  James Orr, LSCE/IPSL, CEA-CNRS-UVSQ, 91191 Gif-sur-Yvette, France
#  15 January 2014

#=======================================================================
#                   define desired precision
#=======================================================================

# set to  2 if you wish results in DOUBLE precision
# set to 1 or 0 if SINGLE
PRECISION = 2
#PRECISION = 1

# mapping between Fortran and C types
ifeq (${PRECISION}, 2)
    KIND_MAP = kind_map_d
else
    KIND_MAP = kind_map_s
endif

#=======================================================================
#=======================================================================

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

FCFLAGS = -fPIC -cpp -DUSE_PRECISION=$(PRECISION)
#DEBUGFLAGS = -g
#LDFLAGS = -L/usr/local/lib -lnetcdf -lnetcdff
LDFLAGS = -L./ -lmocsy
INCLUDEFLAGS = -I/usr/local/include -Isrc


# List of executables to be built within the package
PROGRAMS = libmocsy.a mocsy.so test_mocsy test_errors test_derivauto test_derivnum test_buffesm test_phizero test_kprime test_kzero

# "make" builds all
all: $(PROGRAMS)

#---------------------------------------------------------------------------

# Look for .f90 files first in the 'src' directory, then in the 'examples' directory

vpath %.f90 src
vpath %     examples

#vpath %.h src

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
          vars.f90 \
          derivauto.f90 \
          derivnum.f90 \
          errors.f90 \
          buffesm.f90 \
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
        vars.o \
	derivauto.o \
	derivnum.o \
	errors.o \
        buffesm.o \
        p2fCO2.o \
        f2pCO2.o \
        gasx.o

EXEC = test_mocsy

library = libmocsy.a
#---------------------------------------------------------------------------

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
#%: %.o
#	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

#---------------------------------------------------------------------------
# Build the mocsy library containing the object files (not used, illustration only)
$(library): DNAD.o $(OBJS)
	ar cr $(library) DNAD.o $(OBJS)

# Build the Fortran program executable that tests the mocsy library (test_mocsy)
#<<<<<<< HEAD:Makefile
#$(EXEC): $(EXEC).o $(OBJS) $(library) 
#	$(FC) $(FCFLAGS) -o $@ $@.o $(LDFLAGS)
#
## Build the shared object file for python
#mocsy.so: $(OBJS)
#	cp src/*.f90 .
#	f2py -c $(SOURCES) -m mocsy --fcompiler=gnu95 --f90flags=-O3
#	rm $(SOURCES)
#=======
$(EXEC): $(EXEC).o DNAD.o $(OBJS) test_mocsy.o $(library) 
	$(FC) $(FCFLAGS) -o $@ $@.o $(LDFLAGS)

# Build the shared object file for python
mocsy.so: DNAD.o $(SOURCES)
	cp src/*.f90 .
	# Select the kind map
	cp -f -s src/$(KIND_MAP) .f2py_f2cmap
	f2py -c $(SOURCES) skip: varsolver_dnad : skip: constants_dnad :        \
	skip: sw_ptmp_dnad : skip: sw_temp_dnad : skip: sw_adtg_dnad :          \
	skip: rho_dnad : skip: equation_at_dnad : skip: solve_at_general_dnad : \
	DNAD.o -m mocsy --fcompiler=gnu95 --f90flags="$(FCFLAGS)"
	rm $(SOURCES) DNAD.f90
#---------------------------------------------------------------------------
# Other test programs
#test_errors:  $(LIBSRC_OBJECTS) test_errors.o $(library) 
#	${FC} ${FCFLAGS} -o $@ $@.f90 $(LDFLAGS)

test_errors: test_errors.o $(OBJS) $(library) 
	${FC} ${FCFLAGS} -o $@ $@.o $(LDFLAGS)

test_derivauto: test_derivauto.o $(OBJS) $(library) 
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS)

test_derivnum: test_derivnum.o $(OBJS) $(library) 
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS)

test_buffesm: test_buffesm.o $(OBJS) $(library) 
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS)

test_phizero: test_phizero.o $(OBJS) $(library) 
	${FC} ${FCFLAGS} -o $@ $@.o $(LDFLAGS)

test_kprime: test_kprime.o $(OBJS) $(library) 
	${FC} ${FCFLAGS} -o $@ $@.o $(LDFLAGS)

test_kzero: test_kzero.o $(OBJS) $(library) 
	${FC} ${FCFLAGS} -o $@ $@.o $(LDFLAGS)

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ 

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.so *.a

veryclean: clean
	rm -f *~ $(PROGRAMS) 
