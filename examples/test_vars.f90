!> \file test_mocsy.f90
!! \BRIEF Fortran 90 program to test mocsy.f90
PROGRAM test_vars

   USE msingledouble
   USE mvars

   IMPLICIT NONE

!  For vars routine (called below)
!  "vars" Output variables:
   REAL(kind=rx), DIMENSION(1) :: ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  "vars" Input variables
   INTEGER :: N
   REAL(kind=rx), DIMENSION(1) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat, lon
!  "vars" Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2, optS


!  Local variables:
   INTEGER :: i
   REAL(kind=rx) :: H   ! concentration of [H+]


!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tcsv'    ! input temperature, variable 'temp' is actually conservative temp [°C]
   optP    = 'db'      ! input variable 'depth' is in 'DECIBARS'
   optB    = 'u74'
   optK1K2 = 'l'
   optKf   = 'dg'
   optS    = 'Sabs'
!> Simple input data (with CONCENTRATION units typical for DATA)
!> (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
   DO i = 1,1
     temp(i)   = 18.0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
     sal(i)    = 35.0           !Salinity (practical scale)
     alk(i)    = 2300.*1.e-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
     dic(i)    = 2000.*1.e-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
     sil(i)    = 60.e-6
     phos(i)   = 2.e-6
     depth(i) = 0.
     Patm(i)   = 1.0            !Atmospheric pressure (atm)
     N = i
   END DO

   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                             &  ! INPUT
             optCON='mol/kg', optT='Tcsv   ', optP='db', optB=optB, optK1K2=optK1K2, optKf=optKf, optS='Sabs')


   H = 10**(-ph(1))
   write (*,*) "ph",    ph
   write (*,*) "pco2", pco2
   write (*,*) "fco2", fco2
   write (*,*) "co2",   co2
   write (*,*) "hco3", hco3
   write (*,*) "co3",   co3
   write (*,*) "omegaA", OmegaA
   write (*,*) "omegaC", OmegaC


  STOP
END PROGRAM test_vars

    
