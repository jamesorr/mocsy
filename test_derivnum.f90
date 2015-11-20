!> \file test_derivnum.f90
!! \BRIEF Fortran 90 program to test num_deriv.f90
PROGRAM test_derivnum

   USE msingledouble
   USE mvars
   USE mderivnum

   IMPLICIT NONE


!  Output variables:
   REAL(kind=rx), DIMENSION(1) :: h, ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  derivative of "vars" Output variables:
   REAL(kind=rx), DIMENSION(1) :: dh_dDIC, dpco2_dDIC, dfco2_dDIC, dco2_dDIC, dhco3_dDIC, dco3_dDIC, dOmegaA_dDIC, dOmegaC_dDIC
!  Input variables
   REAL(kind=rx), DIMENSION(1) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
!  Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2

!  Local variables:
   CHARACTER(3) ::  derivar
  !     derivar = 3-character identifier of input variable with respect to which derivative is requested
  !               possibilities are 'alk', 'dic', 'pho', 'sil', 'tem', or 'sal'
  !

!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'm'       ! input variable 'depth' is in meters
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'pf'
!> Simple input data (with CONCENTRATION units typical for DATA)
!> (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
    temp(1)   = 20.0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
    sal(1)    = 35.0           !Salinity (practical scale)
    alk(1)    = 2300.*1.d-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
    dic(1)    = 2000.*1.d-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
    sil(1)    = 60.*1.d-6
    phos(1)   = 2.*1.d-6
    depth(1)  = 0.
    Patm(1)   = 1.0            !Atmospheric pressure (atm)
    lat(1)    = 0.

!  Select input var 'x', choosing set of dy_i/dx to be computed, where y_i are the diff output vars
   derivar = 'dic'
!  derivar = 'alk'
!  derivar = 'pho'
!  derivar = 'sil'
!  derivar = 'tem'
!  derivar = 'sal'
   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                             &  ! INPUT
             optCON, optT, optP    )
   h(1) = 10**(- ph(1))

   call derivnum (dh_dDIC, dpco2_dDIC, dfco2_dDIC, dco2_dDIC, dhco3_dDIC,             &
                  dco3_dDIC, dOmegaA_dDIC, dOmegaC_dDIC,                              &
                  temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1, derivar,             &
                  optCON, optT, optP                                                )

   write (*,*) "Absolute derivatives"
   
   write (*,*) "dh_dDIC",      dh_dDIC(1)
   write (*,*) "dpco2_dDIC",   dpco2_dDIC(1)
   write (*,*) "dfco2_dDIC",   dfco2_dDIC(1)
   write (*,*) "dco2_dDIC",    dco2_dDIC(1)
   write (*,*) "dhco3_dDIC",   dhco3_dDIC(1)
   write (*,*) "dco3_dDIC",    dco3_dDIC(1)
   write (*,*) "dOmegaA_dDIC", dOmegaA_dDIC(1)
   write (*,*) "dOmegaC_dDIC", dOmegaC_dDIC(1)

   write (*,*) "Relative derivatives"
   
   write (*,*) "dlogh_dlogDIC",      dh_dDIC(1) * dic / h
   write (*,*) "dlogpco2_dlogDIC",   dpco2_dDIC(1) * dic / pco2
   write (*,*) "dlogfco2_dlogDIC",   dfco2_dDIC(1) * dic / fco2
   write (*,*) "dlogco2_dlogDIC",    dco2_dDIC(1) * dic / co2
   write (*,*) "dloghco3_dlogDIC",   dhco3_dDIC(1) * dic / hco3
   write (*,*) "dlogco3_dlogDIC",    dco3_dDIC(1) * dic / co3
   write (*,*) "dlogOmegaA_dlogDIC", dOmegaA_dDIC(1) * dic / OmegaA
   write (*,*) "dlogOmegaC_dlogDIC", dOmegaC_dDIC(1) * dic / OmegaC

  STOP
END PROGRAM test_derivnum

    
