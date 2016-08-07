!> \file test_mocsy.f90
!! \BRIEF Fortran 90 program to test mocsy.f90
PROGRAM test_mocsy

   USE msingledouble
   USE mconstants
   USE mvars
   USE mderivauto

   IMPLICIT NONE

!  For vars routine (called below)
!  "vars" Output variables:
   REAL(kind=rx), DIMENSION(100) :: ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  "vars" Input variables
   INTEGER :: N
   REAL(kind=rx), DIMENSION(100) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
   REAL(kind=rx), DIMENSION(6,100) :: ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, hco3_deriv, co3_deriv, OmegaA_deriv, OmegaC_deriv
   REAL(kind=rx) ::  gamma_DIC, gamma_Alk, beta_DIC, beta_Alk, omega_DIC, omega_Alk
!  "vars" Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2

  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  REAL(kind=r8), DIMENSION(6) :: K0
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), DIMENSION(6) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), DIMENSION(6) :: K2
  !> equilibrium constant for dissociation of boric acid 
  REAL(kind=r8), DIMENSION(6) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  REAL(kind=r8), DIMENSION(6) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  REAL(kind=r8), DIMENSION(6) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride 
  !! either from Dickson and Riley (1979) or from Perez and Fraga (1987), depending on optKf
  REAL(kind=r8), DIMENSION(6) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  REAL(kind=r8), DIMENSION(6) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  REAL(kind=r8), DIMENSION(6) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), DIMENSION(6) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), DIMENSION(6) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), DIMENSION(6) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  REAL(kind=r8), DIMENSION(6) :: Ksi
  !> total sulfate (Morris & Riley, 1966)
  REAL(kind=r8), DIMENSION(6) :: St
  !> total fluoride  (Riley, 1965)
  REAL(kind=r8), DIMENSION(6) :: Ft
  !> total boron
  !! from either Uppstrom (1974) or Lee et al. (2010), depending on optB
  REAL(kind=r8), DIMENSION(6) :: Bt

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
             optCON='mol/kg', optT='Tinsitu', optP='db', optB='l10', optK1K2=optK1K2, optKf='dg')     ! INPUT OPTIONS

   call derivauto(ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, hco3_deriv, co3_deriv,         &
             OmegaA_deriv, OmegaC_deriv,                                                      &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                             &  ! INPUT
             optCON='mol/kg', optT='Tinsitu', optP='db', optB='l10', optK1K2=optK1K2, optKf='dg')     ! INPUT OPTIONS
             
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
   write(*,51)

!  Need to change to call latest version of constants & derivnum (to get constand 
!  call constants (K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa,                   &
!                     K1p, K2p, K3p, Ksi,                                    &
!                     St, Ft, Bt,                                            &
!                     temp, sal, Patm,                                       &
!                     depth, lat, 6,                                         &
!                     optT='Tinsitu', optP='db', optB='l10', optK1K2=optK1K2, optKf='dg',  &
!                     K0_deriv=K0_deriv, Kb_deriv=Kb_deriv, Kspa_deriv=Kspa_deriv )
!   
!!  Print out derivatives of thermodynamic constants
!   write(*,50)
!   write(*,300)
!   do i=1,N
!     write(*,301) K0(i), K0_deriv(1,i), K0_deriv(2,i), Kb(i), Kb_deriv(1,i), Kb_deriv(2,i),    &
!        Kspa(i), Kspa_deriv(1,i), Kspa_deriv(2,i)
!   end do
!   write(*,52)
   
  50 format(/)
  51 format(179('-'))
! 52 format(106('-'))


 200 format('Typical DATA output',                                                                              / &
            179('-'),                                                                                                    / &
            '                                                                               in situ         in situ',    / &
            '  pH     pCO2   fCO2     CO2*       HCO3-       CO32-      OmegaA OmegaC  R    Density Press  Temperature', &
            ' gamma_DIC  gamma_Alk     beta_DIC     beta_Alk    omega_DIC     omega_Alk', / &
            '(total) (uatm) (uatm)  (mol/kg)    (mol/kg)     (mol/kg)                       (kg/m3) (db)      (C)',      / &
            179('-')                                                                                                     )

 201 format(f7.4, 2f7.1, 3(e12.4), 3f7.2, f8.2, f7.1, f8.3, 6f13.9)

! 300 format('Constants and their derivatives',                                                   / &
!            106('-'),                                                                                                    / &
!            '  K0        dK0/dT       dK0/dS       Kb        dKb/dT      dKb/dS       Kspa       dKspa/dT    dKspa/dT', / &
!            106('-')                                                                                                     )
!
! 301 format(f10.7, 8(e12.4, e12.4) )

  STOP
END PROGRAM test_mocsy

