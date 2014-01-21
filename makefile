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
#INCLUDEFLAGS = -I/usr/local/include

#---------------------------------------------------------------------------
MOCSY = singledouble.f90 \
       sw_adtg.f90 \
       sw_ptmp.f90 \
       sw_temp.f90 \
       tpot.f90 \
       tis.f90 \
       p80.f90 \
       rho.f90 \
       rhoinsitu.f90 \
       depth2press.f90 \
       constants.f90 \
       vars.f90 

OBJS = $(MOCSY) \
       test_mocsy.f90 

EXEC = test_mocsy
#---------------------------------------------------------------------------

# use Pattern rules (see gnumake documentation)
%.o: %.f
	$(FC) $(FCFLAGS) -c -o $@ $(INCLUDEFLAGS) $*.f

%.o: %.f90
	$(F90) $(FCFLAGS) -c -o $@ $(INCLUDEFLAGS) $*.f90

#---------------------------------------------------------------------------
$(EXEC): $(OBJS) 
	$(FC) $(FCFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS) $(INCLUDEFLAGS)

$(EXEC).o: makefile
#---------------------------------------------------------------------------

clean:  
	rm *.mod *.so
	rm $(EXEC) 
	#rm $(EXEC) $(EXEC2)

python: 
	f2py -c $(MOCSY) -m mocsy --fcompiler=gnu95 --f90flags=-O3
	
