!> \file test_phizero.f90
!! \BRIEF Fortran 90 program to test phizero routine in gasx.f90
PROGRAM test_phizero

   USE msingledouble
   USE gasx

   IMPLICIT NONE

   INTEGER, PARAMETER :: n = 1
   
!  Output variables:
   REAL(kind=r8), DIMENSION(1) :: phi0_cfc11, phi0_cfc12, phi0_sf6, phi0_co2, phi0_n2o
   
!  Input variables
   REAL(kind=rx), DIMENSION(1) :: temp, sal
   
!  Input at standard T and S
!  temp(1)   = 25.0    
!  sal(1)    = 35.0    

   temp(1)   = 4.0    
   sal(1)    = 34.0    

!  Call
   call phizero('cfc11', temp, sal, n, phi0_cfc11)
   call phizero('cfc12', temp, sal, n, phi0_cfc12)
   call phizero('sf6',   temp, sal, n, phi0_sf6)
   call phizero('co2',   temp, sal, n, phi0_co2)
   call phizero('n2o',   temp, sal, n, phi0_n2o)

   write (*,*) "phi0_cfc11 = ", phi0_cfc11
   write (*,*) "phi0_cfc12 = ", phi0_cfc12
   write (*,*) "phi0_sf6   = ", phi0_sf6
   write (*,*) "phi0_co2   = ", phi0_co2
   write (*,*) "phi0_n2o   = ", phi0_n2o

   write (*, *) ' '
   write (*,*) 'Results above should be 2.136063e-2, 5.370289e-3, 3.560764e-4, 5.517209e-2, 4.080887e-2, respectively:'
   write (*,*) '-> these values are essentially the same as in the papers from Warner & Weiss (1985), Bullister et al. (2002), and' 
   write (*,*) '   Weiss & Price (1980), after correcting for density (using 1.028).'

   STOP
END PROGRAM test_phizero

    
