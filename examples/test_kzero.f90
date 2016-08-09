!> \file test_kzero.f90
!! \BRIEF Fortran 90 program to test kzero routine in gasx.f90
PROGRAM test_kzero

   USE msingledouble
   USE gasx

   IMPLICIT NONE

   INTEGER, PARAMETER :: n = 1
   
!  Computed variables:
   REAL(kind=r8), DIMENSION(1) :: k0_co2, k0_n2o
   
!  Input variables
   REAL(kind=rx), DIMENSION(1) :: temp, sal
   
!  Input at standard T and S
!  temp(1)   = 25.0    
!  sal(1)    = 35.0    

   temp(1)   = 4.0    
   sal(1)    = 34.0    

   call kzero('co2',   temp, sal, n, k0_co2)
   call kzero('n2o',   temp, sal, n, k0_n2o)

   write (*,*) "k0_co2   = ", k0_co2
   write (*,*) "k0_n2o   = ", k0_n2o

   write (*, *) ' '
   write (*,*) 'Results above should be 5.585e-2 and 4.132e-2 respectively:'
   write (*,*) '-> these values are the same as those in from Weiss (1974, Table 2) and' 
   write (*,*) '   Weiss & Price (1980, Table 3), respectively, at T=4 and S=34'

   STOP
END PROGRAM test_kzero

    
