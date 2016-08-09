!> \file test_phizero.f90
!! \BRIEF Fortran 90 program to test phizero routine in gasx.f90
PROGRAM test_kprime

   USE msingledouble
   USE gasx

   IMPLICIT NONE

   INTEGER, PARAMETER :: n = 1
   
!  Output variables:
   REAL(kind=r8), DIMENSION(1) :: kp_cfc11, kp_cfc12, kp_sf6
   
!  Input variables
   REAL(kind=rx), DIMENSION(1) :: temp, sal
   
!  Input at standard T and S
!  temp(1)   = 25.0    
!  sal(1)    = 35.0    

   temp(1)   = 4.0    
   sal(1)    = 34.0    

!  Call
   call kprime('cfc11', temp, sal, n, kp_cfc11)
   call kprime('cfc12', temp, sal, n, kp_cfc12)
   call kprime('sf6',   temp, sal, n, kp_sf6)

   write (*,*) "kp_cfc11 = ", kp_cfc11
   write (*,*) "kp_cfc12 = ", kp_cfc12
   write (*,*) "kp_sf6   = ", kp_sf6

   write (*, *) ' '
   write (*,*) 'Results above should be very close to 2.156e-2, 5.4186e-3, and 3.5842e-4, respectively:'
   write (*,*) '-> these values are about the same as those in from Warner & Weiss (1985, Tables 3 & 4), and'
   write (*,*) '   Bullister et al. (2002, Table 2), at T=4 and S=34.'
   write (*,*) '   These first two check values are those from Warner & Weiss in mol kg-1 atm-1 multiplied by 1.028.'
   write (*,*) '   These 3rd one is identical to that from Bullister et al., but their units given are wrong (Table 2 footnote).'


   STOP
END PROGRAM test_kprime

    
