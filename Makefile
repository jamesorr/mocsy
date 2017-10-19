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
#               define GSW-fortran files 
#=======================================================================
GSW = src/GSW

GSW_MOD_SRCS := \
	$(GSW)/gsw_mod_kinds.f90 \
	$(GSW)/gsw_mod_teos10_constants.f90 \
	$(GSW)/gsw_mod_toolbox.f90 \
	$(GSW)/gsw_mod_error_functions.f90 \
	$(GSW)/gsw_mod_baltic_data.f90 \
	$(GSW)/gsw_mod_saar_data.f90 \
	$(GSW)/gsw_mod_specvol_coefficients.f90

GSW_MOD_OBJS := $(GSW_MOD_SRCS:.f90=.o)

GSW_TOOL_SRCS := \
	$(GSW)/gsw_t_from_ct.f90 \
	$(GSW)/gsw_ct_from_t.f90 \
	$(GSW)/gsw_ct_from_pt.f90 \
	$(GSW)/gsw_pt_from_ct.f90 \
	$(GSW)/gsw_pt_from_t.f90 \
	$(GSW)/gsw_pt0_from_t.f90 \
	$(GSW)/gsw_gibbs_pt0_pt0.f90 \
	$(GSW)/gsw_entropy_part.f90 \
	$(GSW)/gsw_entropy_part_zerop.f90 \
	$(GSW)/gsw_gibbs.f90 \
	$(GSW)/gsw_sp_from_sa.f90 \
	$(GSW)/gsw_sa_from_sp.f90 \
	$(GSW)/gsw_saar.f90 \
	$(GSW)/gsw_sp_from_sa_baltic.f90 \
	$(GSW)/gsw_sa_from_sp_baltic.f90 \
	$(GSW)/gsw_util_xinterp1.f90 \
	$(GSW)/gsw_util_indx.f90 \
	$(GSW)/gsw_add_barrier.f90 \
	$(GSW)/gsw_add_mean.f90 \
	$(GSW)/gsw_rho.f90 \
	$(GSW)/gsw_specvol.f90

GSW_TOOL_OBJS := $(GSW_TOOL_SRCS:.f90=.o)

#=======================================================================
#=======================================================================

# To use another Fortran compiler, replace "f95" in FC and F90 lines with your compiler command
# For example, comment out 2 lines below with "f95" & uncomment the following 2 lines for "gfortran"
#FC = fort77
#FC = xlf
#FC = f95
#F90 = f95
FC = gfortran -ffree-line-length-none
F90 = gfortran -ffree-line-length-none
#FC = ifort
#F90 = ifort

FCFLAGS = -fPIC -cpp -DUSE_PRECISION=$(PRECISION)
#DEBUGFLAGS = -g
LDFLAGS = -L./ -lmocsy
INCLUDEFLAGS = -Isrc


# List of executables to be built within the package
PROGRAMS = libmocsy.a mocsy.so test_mocsy test_errors test_derivauto test_derivnum test_buffesm test_phizero test_kprime test_kzero

# "make" builds all
all: $(PROGRAMS)

#---------------------------------------------------------------------------

# Look for .f90 test files in the 'examples' directory
vpath %     examples

# Attention: src/singledouble.f90 is automatically generated
src/singledouble.f90 : src/singledouble.m4
	m4 -DUSE_PRECISION=$(PRECISION) $^ > $@

SOURCES = src/singledouble.f90 \
          src/eos.f90 \
          src/sw_adtg.f90 \
          src/sw_ptmp.f90 \
          src/sw_temp.f90 \
          src/tpot.f90 \
          src/tis.f90 \
          src/p80.f90 \
          src/phsolvers.f90 \
          src/rho.f90 \
          src/rhoinsitu.f90 \
          src/depth2press.f90 \
          src/constants.f90 \
          src/varsolver.f90 \
          src/vars.f90 \
          src/derivauto.f90 \
          src/derivnum.f90 \
          src/errors.f90 \
          src/buffesm.f90 \
          src/p2fCO2.f90 \
          src/f2pCO2.f90 \
          src/gasx.f90 

OBJS := $(SOURCES:.f90=.o)

EXEC = test_mocsy

library = libmocsy.a
#---------------------------------------------------------------------------

# Build the mocsy library containing the object files
$(library): $(GSW_MOD_OBJS) $(GSW_TOOL_OBJS) src/DNAD.o $(OBJS)
	ar cr $@ $^

# Build the Fortran program executable that tests the mocsy library (test_mocsy)
$(EXEC): $(library) $(EXEC).o
	$(FC) $(FCFLAGS) -o $@ $@.o $(LDFLAGS)

# Build the shared object file for python
mocsy.so: $(GSW_MOD_OBJS) $(GSW_TOOL_OBJS) src/DNAD.o $(SOURCES)
	# cp src/*.f90 .
	# Select the kind map
	cp -f -s src/$(KIND_MAP) .f2py_f2cmap
	f2py -c -L. $(SOURCES) skip: vars_sprac : skip: vars_pertK : skip: varsolver_dnad :   \
	    skip: constants_dnad : skip: sw_ptmp_dnad : skip: sw_temp_dnad : skip: sw_adtg_dnad :   \
	    skip: rho_dnad : skip: equation_at_dnad : skip: solve_at_general_dnad :                 \
	    $(GSW_MOD_OBJS) $(GSW_TOOL_OBJS) src/DNAD.o -m mocsy --fcompiler=gnu95                      \
	    --f90flags="$(FCFLAGS) $(INCLUDEFLAGS)"
	# rm $(SOURCES) DNAD.f90
#---------------------------------------------------------------------------
# Other test programs

test_vars: $(library) test_vars.o
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS) 

test_errors: $(library) test_errors.o
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS) 

test_derivauto: $(library) test_derivauto.o
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS)

test_derivnum: $(library) test_derivnum.o
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS)

test_buffesm: $(library) test_buffesm.o
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS)

test_phizero: $(library) test_phizero.o
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS) 

test_kprime: $(library) test_kprime.o
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS) 

test_kzero: $(library) test_kzero.o
	${F90} ${FCFLAGS} -o $@ $@.o $(LDFLAGS) 

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ 

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

%.o: %.F90
	$(FC) $(FCFLAGS) $(INCLUDEFLAGS) -c $< -o $@

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o src/*.o $(GSW)/*.o *.mod *.so *.a src/singledouble.f90

veryclean: clean
	rm -f *~ $(PROGRAMS) 
