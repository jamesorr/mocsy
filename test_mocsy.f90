!> \file test_mocsy.f90
!! \BRIEF Fortran 90 program to test mocsy.f90
PROGRAM test_mocsy

   USE msingledouble
   USE mvars

   IMPLICIT NONE


!  For vars routine (called below)
!  "vars" Output variables:
   REAL(kind=r4), DIMENSION(100) :: ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  "vars" Input variables
   INTEGER :: N
   REAL(kind=r4), DIMENSION(100) :: tempot, sal, alk, dic, sil, phos, depth, lat
!  "vars" Input options
   CHARACTER(10) :: optRHO, optT, optP, optB, optKf, optK1K2

!  Local variables:
   INTEGER :: i

!> Simple input data (with CONCENTRATION units typical for MODELS)
   DO i = 1,6
     tempot(i) = 2.0
     sal(i)    = 35.0
     alk(i)    = 2295.*1.025e-3 !Rough convert of S. Ocean ave surf ALK (umol/kg) to model units (mol/m3)
     dic(i)    = 2154.*1.025e-3 !Rough convert of S. Ocean ave surf DIC (umol/kg) to model units (mol/m3)
     sil(i)    = 50.  *1.0e-3   !Convert observed S. Ocean ave surf SiO2 (umol/L) to model units (mol/m3) 
     phos(i)   = 1.8  *1.0e-3   !Convert observed S. Ocean ave surf PO4  (umol/L) to model units (mol/m3)
     lat(i)    = -60.           !Latitude used to convert model depth to pressure (Saunders, 1981, JPO)
     depth(i) = real(i-1) * 1000. ! depth varies from 0 to 5000 m by 1000 m
     N = i
   END DO

!> OPTIONS: see complete documentation in 'vars' subroutine
!> Typical options for MODELS
   optRHO  = 'mol/m3'   ! input concentrations are in MOL/M3
   optT    = 'Tpot'     ! input temperature, variable 'tempot', is POTENTIAL temp [°C]
   optP    = 'm'        ! input variable 'depth' is in METERS
   optB    = 'l10'      ! Lee et al. (2010) formulation for total boron
   optK1K2 = 'l'          ! Lueker et al. (2000) formulations for K1 & K2 (best practices)
   optKf   = 'dg'       ! Dickson & Riley (1979) formulation for Kf (recommeded by Dickson & Goyet, 1994)

!> Call mocsy's main subroutine to compute carbonate system variables: pH, pCO2, fCO2, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R
!> FROM temperature, salinity, total alkalinity, dissolved inorganic carbon, silica, phosphate, depth (or pressure) (1-D arrays)
   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &  ! OUTPUT
             tempot, sal, alk, dic, sil, phos, depth, lat, N,                          &  ! INPUT
             optRHO, optT, optP, optB, optK1K2, optKf)                                    ! INPUT OPTIONS
!  Print out results
   write(*,50)
   write(*,100)
   do i=1,N
     write(*,101) ph(i), pco2(i), fco2(i), co2(i), hco3(i), co3(i),  OmegaA(i), OmegaC(i), BetaD(i), rhoSW(i), p(i), tempis(i)
   end do
   write(*,51)

!  write(*,*) (dic(i), i=1,6)

!> Repeat the operation, but with options for observations
!> Typical options for observations
   optRHO  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'tempot' is actually IN SITU temp [°C]
   optP    = 'db'      ! input variable 'depth' is in 'DECIBARS'
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'dg'
!> Simple input data (with CONCENTRATION units typical for DATA)
!> (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
   DO i = 1,6
     alk(i)    = 2295.*1.e-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
     dic(i)    = 2154.*1.e-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
     sil(i)    = 50.  /1.025e6    ! Rough convert obs. S. Ocean ave surf SiO2 (umol/L) to mocsy data units (mol/kg) 
     phos(i)   = 1.8  /1.025e6    ! Rough convert obs. S. Ocean ave surf PO4  (umol/L) to model data units (mol/kg)
     depth(i) = real(i-1) * 1000. ! Vary depth from 0 to 5000 db by 1000 db
     N = i
   END DO

   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
             tempot, sal, alk, dic, sil, phos, depth, lat, N,                                 &  ! INPUT
             optRHO='mol/kg', optT='Tinsitu', optP='db', optB='l10', optK1K2='m10', optKf='dg')  ! INPUT OPTIONS
!  Print out results (typical for data: concentration units differ)
   write(*,50)
   write(*,200)
   do i=1,N
     write(*,201) ph(i), pco2(i), fco2(i), co2(i), hco3(i), co3(i),  OmegaA(i), OmegaC(i), BetaD(i), rhoSW(i), p(i), tempis(i)
   end do
   write(*,51)
   write(*,202)
   write(*,203)
   write(*,50)

  50 format(/)
  51 format(105('-'))

 100 format('Table 1: Typical MODEL output',                                                                             / &
            105('-'),                                                                                                    / &
            '                                                                               in situ         in situ',    / &
            '  pH     pCO2   fCO2     CO2*       HCO3-       CO32-      OmegaA OmegaC  R    Density Press  Temperature', / &
            '(total) (uatm) (uatm)  (mol/m3)    (mol/m3)    (mol/m3)                        (kg/m3) (db)      (C)',      / &
            105('-')                                                                                                     )

 200 format('Table 2: Typical DATA output',                                                                              / &
            105('-'),                                                                                                    / &
            '                                                                               in situ         in situ',    / &
            '  pH     pCO2   fCO2     CO2*       HCO3-       CO32-      OmegaA OmegaC  R    Density Press  Temperature', / &
            '(total) (uatm) (uatm)  (mol/kg)    (mol/kg)     (mol/kg)                       (kg/m3) (db)      (C)',      / &
            105('-')                                                                                                     )

 101 format(f7.4, 2f7.1, 3(e12.4), 3f7.2, f8.2, f7.1, f8.3)
 201 format(f7.4, 2f7.1, 3(e12.4), 3f7.2, f8.2, f7.1, f8.3)
 202 format('*Table 2 differs slightly from Table 1 because input DIC, Alk, pressure, and temperature differ.')
 203 format('*Concentration units differ between Tables 1 and 2 because input options differ')

  STOP
END PROGRAM test_mocsy

