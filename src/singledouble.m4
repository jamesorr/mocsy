dnl 
dnl This file is a template file for single_double.f90 : it must be processed by m4
dnl -------------------------------------------------------------------------------

!> \file singledouble.f90
!! \BRIEF 
!> Module that defines single and double precision - used by all other modules
MODULE msingledouble
! INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)
! INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: rx = KIND(ifelse(USE_PRECISION,2,1.0d0,1.0e0))

  INTEGER, PARAMETER :: r8 = KIND(1.0d0)
  INTEGER, PARAMETER :: wp = KIND(1.0d0)
END MODULE msingledouble