!> \file test_mocsy.f90
!! \BRIEF Fortran 90 program to test mocsy.f90
PROGRAM test_derivauto

   USE msingledouble
   USE mderivauto
   USE mvars

   IMPLICIT NONE


!  For vars routine (called below)
!  "vars" Output variables:
   REAL(kind=rx), DIMENSION(1) :: ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  "vars" Input variables
   INTEGER :: N
   REAL(kind=rx), DIMENSION(1) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
   REAL(kind=rx), DIMENSION(6,1) :: ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, &
                     hco3_deriv, co3_deriv, omegaa_deriv, omegac_deriv
!  "vars" Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2


!  Local variables:
   INTEGER :: i
   REAL(kind=r8) :: H
   REAL(kind=r8), DIMENSION(6) :: H_deriv
   CHARACTER*4 :: invar(6)
   
   invar(1) = 'Alk '
   invar(2) = 'DIC '
   invar(3) = 'Phos'
   invar(4) = 'Sil '
   invar(5) = 'T   '
   invar(6) = 'S   '

!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'm'      ! input variable 'depth' is in 'DECIBARS'
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'dg'
!> Simple input data (with CONCENTRATION units typical for DATA)
!> (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
   DO i = 1,1
     temp(i)   = 18.0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
     sal(i)    = 35.0           !Salinity (practical scale)
     alk(i)    = 2300.*1.e-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
     dic(i)    = 2000.*1.e-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
     sil(i)    = 60.*1.e-6
     phos(i)   = 2.*1.e-6
     depth(i)  = 0.
     Patm(i)   = 1.0            !Atmospheric pressure (atm)
     lat(i)    = 0.
     N = i
   END DO

   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                             &  ! INPUT
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf    )

   call derivauto(ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, hco3_deriv, co3_deriv,   &
                omegaa_deriv, omegac_deriv,                                             &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                    &
                optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf )

    ! [H+] concentration
    H = 10.0**(-ph(1))
    ! derivative of [H+] deduced from that of pH
    H_deriv = - H * ph_deriv(:,1) * log(10.0)
    
    write (*,*) "derivatives :"
    write (*,*) "===========  "

    write (*,*) "      dH/dx         dpCO2/dx       dfCO2/dx", &
          "       d[CO2*]/dx      d[HCO3-]/dx    d[CO3--]/dx     dOmegaA/dx     dOmegaC/dx"
    
    DO i = 1,6
       write (*,"(A3, 8ES15.6)") invar(i), H_deriv(i), pco2_deriv(i,1), fco2_deriv(i,1), co2_deriv(i,1), hco3_deriv(i,1), &
       co3_deriv(i,1), omegaa_deriv(i,1), omegac_deriv(i,1)
    END DO
    

  STOP
END PROGRAM test_derivauto

