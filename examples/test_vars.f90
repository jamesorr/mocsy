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
   REAL(kind=rx), DIMENSION(1) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
!  REAL(kind=rx), DIMENSION(6,1) :: ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, &
!                    hco3_deriv, co3_deriv, omegaa_deriv, omegac_deriv
!  "vars" Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2


!  Local variables:
   INTEGER :: i
   REAL(kind=rx) :: H   ! concentration of [H+]


!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'db'      ! input variable 'depth' is in 'DECIBARS'
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'dg'
!> Simple input data (with CONCENTRATION units typical for DATA)
!> (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
   DO i = 1,1
     temp(i)   = 20.0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
     sal(i)    = 35.0           !Salinity (practical scale)
     alk(i)    = 2300.*1.e-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
     dic(i)    = 2000.*1.e-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
     sil(i)    = 0.
     phos(i)   = 0.
     depth(i) = 0.
     Patm(i)   = 1.0            !Atmospheric pressure (atm)
     N = i
   END DO

!  call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
!            temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                             &  ! INPUT
!            optCON='mol/kg', optT='Tinsitu', optP='db', &
!            pco2_deriv=pco2_deriv)

!  write (*,*) "pco2_deriv[1,0]", pco2_deriv(2,1)
   
!   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
!             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                             &  ! INPUT
!             optCON='mol/kg', optT='Tinsitu', optP='db', &
!             ph_deriv=ph_deriv, pco2_deriv=pco2_deriv, fco2_deriv=fco2_deriv, co2_deriv=co2_deriv, &
!             hco3_deriv=hco3_deriv, co3_deriv=co3_deriv, omegaa_deriv=omegaa_deriv, omegac_deriv=omegac_deriv)

   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                             &  ! INPUT
             optCON='mol/kg', optT='Tinsitu', optP='db')


   write (*,*) "pco2_deriv[1,0]", pco2_deriv(2,1)
   
   H = 10**(-ph(1))
   write (*,*) "h_deriv",    -H * ph_deriv * log(10.0)
   write (*,*) "pco2_deriv", pco2_deriv
   write (*,*) "fco2_deriv", fco2_deriv
   write (*,*) "co2_deriv",   co2_deriv
   write (*,*) "hco3_deriv", hco3_deriv
   write (*,*) "co3_deriv",   co3_deriv
   write (*,*) "omegaA_deriv", omegaa_deriv
   write (*,*) "omegaC_deriv", omegac_deriv


  STOP
END PROGRAM test_vars

    
