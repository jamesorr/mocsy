!> \file test_mocsy.f90
!! \BRIEF Fortran 90 program to test mocsy.f90
PROGRAM test_mocsy

   USE msingledouble
   USE mvars

   IMPLICIT NONE


!  For vars routine (called below)
!  "vars" Output variables:
   REAL(kind=rx), DIMENSION(100) :: ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  "vars" Input variables
   INTEGER :: N
   REAL(kind=rx), DIMENSION(100) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
   REAL(kind=rx), DIMENSION(4,100) :: ph_deriv, pco2_deriv, OmegaA_deriv
   REAL(kind=rx) ::  gamma_DIC, gamma_Alk, beta_DIC, beta_Alk, omega_DIC, omega_Alk
!  "vars" Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2

!  Local variables:
   INTEGER :: i


!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'db'      ! input variable 'depth' is in 'DECIBARS'
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'dg'
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
             optCON='mol/kg', optT='Tinsitu', optP='db', optB='l10', optK1K2=optK1K2, optKf='dg',&
             ph_deriv=ph_deriv, pco2_deriv=pco2_deriv, OmegaA_deriv=OmegaA_deriv)     ! INPUT OPTIONS
!  Print out results (typical for data: concentration units differ)
   write(*,50)
   write(*,200)
   do i=1,N
     gamma_DIC = pco2(i) / pco2_deriv(2,i)
     gamma_Alk = pco2(i) / pco2_deriv(1,i)

     beta_DIC  = -1. / (LOG(10.) * ph_deriv(2,i))
     beta_Alk  = -1. / (LOG(10.) * ph_deriv(1,i))
     
     ! Here, we use Omega of Aragonite (use of Calcite would have been equaly valid)
     omega_DIC = OmegaA(i) / OmegaA_deriv(2,i)
     omega_Alk = OmegaA(i) / OmegaA_deriv(1,i)
     
     write(*,201) ph(i), pco2(i), fco2(i), co2(i), hco3(i), co3(i), &
     OmegaA(i), OmegaC(i), BetaD(i), rhoSW(i), p(i), tempis(i),     &
     gamma_DIC, gamma_Alk, beta_DIC, beta_Alk, omega_DIC, omega_Alk
   end do

  50 format(/)
  51 format(105('-'))


 200 format('Typical DATA output',                                                                              / &
            105('-'),                                                                                                    / &
            '                                                                               in situ         in situ',    / &
            '  pH     pCO2   fCO2     CO2*       HCO3-       CO32-      OmegaA OmegaC  R    Density Press  Temperature', &
            ' gamma_DIC  gamma_Alk  beta_DIC  beta_Alk  omega_DIC omega_Alk', / &
            '(total) (uatm) (uatm)  (mol/kg)    (mol/kg)     (mol/kg)                       (kg/m3) (db)      (C)',      / &
            105('-')                                                                                                     )

 201 format(f7.4, 2f7.1, 3(e12.4), 3f7.2, f8.2, f7.1, f8.3, 6f13.9)

  STOP
END PROGRAM test_mocsy

