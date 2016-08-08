!> \file test_mocsy.f90
!! \BRIEF Fortran 90 program to test mocsy.f90
PROGRAM test_buffesm

   USE msingledouble
   USE mvars
   USE mbuffesm
   
   IMPLICIT NONE


!  For vars routine (called below)
!  "vars" Output variables:g
   REAL(kind=rx), DIMENSION(100) :: ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  "vars" Input variables
   INTEGER :: N
   REAL(kind=rx), DIMENSION(100) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
   REAL(kind=rx), DIMENSION(100)::  gammaDIC, gammaAlk, betaDIC, betaAlk, omegaDIC, omegaAlk, Rf
!  "vars" Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2, optGAS
!  CHARACTER(7) :: optGAS

!  Local variables:
   INTEGER :: i


!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'db'      ! input variable 'depth' is in 'DECIBARS'
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'dg'
   optGAS  = 'Pinsitu'
!> Simple input data (with CONCENTRATION units typical for DATA)
!> (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
   DO i = 1,6
     temp(i)   = 2.0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
     sal(i)    = 35.0           !Salinity (practical scale)
     alk(i)    = 2295.*1.e-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
     dic(i)    = 2154.*1.e-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
     sil(i)    = 0.
     phos(i)   = 0.
     depth(i) = real(i-1) * 1000. ! Vary depth from 0 to 5000 db by 1000 db
     Patm(i)   = 1.0            !Atmospheric pressure (atm)
     N = i
   END DO

   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                             &  ! INPUT
             optCON='mol/kg', optT='Tinsitu', optP='db', optB='l10', optK1K2=optK1K2,         &  ! OPTIONS
             optKf='dg', optGAS=optGAS)                                                          

   call buffesm(gammaDIC, betaDIC, omegaDIC, gammaALK, betaALK, omegaALK, Rf,                    &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                                 &  ! INPUT
             optCON='mol/kg', optT='Tinsitu', optP='db', optB='l10', optK1K2=optK1K2,             &  ! OPTIONS
             optKf='dg', optGAS=optGAS)                                                              

!  Print out results (typical for data: concentration units differ)
   write(*,50)
   write(*,200)
   do i=1,N
     write(*,201) ph(i), pco2(i), fco2(i), co2(i), hco3(i), co3(i), &
     OmegaA(i), OmegaC(i), BetaD(i), rhoSW(i), p(i), tempis(i),     & 
     gammaDIC(i), gammaAlk(i), betaDIC(i), betaAlk(i), omegaDIC(i), omegaAlk(i)
   end do
   write(*,51)

  50 format(/)
  51 format(179('-'))
  52 format(106('-'))


 200 format('Typical DATA output',                                                                              / &
            179('-'),                                                                                                    / &
            '                                                                               in situ         in situ',    / &
            '  pH     pCO2   fCO2     CO2*       HCO3-       CO32-      OmegaA OmegaC  R    Density Press  Temperature', &
            ' gamma_DIC  gamma_Alk     beta_DIC     beta_Alk     omega_DIC    omega_Alk', / &
            '(total) (uatm) (uatm)  (mol/kg)    (mol/kg)     (mol/kg)                       (kg/m3) (db)      (C)',      / &
            179('-')                                                                                                     )

 201 format(f7.4, 2f7.1, 3(e12.4), 3f7.2, f8.2, f7.1, f8.3, 6f13.9)

  STOP
END PROGRAM test_buffesm

