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
   REAL(kind=rx), DIMENSION(1) :: dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx, dco3_dx, dOmegaA_dx, dOmegaC_dx
!  Input variables
   REAL(kind=rx), DIMENSION(1) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
!  Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2

!  Local variables:
   CHARACTER*3, DIMENSION(13) ::  devar = (/'alk','dic','pho','sil','tem','sal','k0 ','k1 ','k2 ','kb ','kw ','ka ','kc '/)
   INTEGER ::  i
   
  !     derivar = 3-character identifier of input variable with respect to which derivative is requested
  !               possibilities are 'alk', 'dic', 'pho', 'sil', 'tem', or 'sal'
  !

!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'db'       ! input variable 'depth' is in meters
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'dg'
!> Simple input data (with CONCENTRATION units typical for DATA)
!> (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
    temp(1)   = 18.0d0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
    sal(1)    = 35.0d0           !Salinity (practical scale)
    alk(1)    = 2300.0d-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
    dic(1)    = 2000.0d-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
    sil(1)    = 0.0d0   ! 60.d-06
    phos(1)   = 0.0d0   !  2.d-06
    sil(1)    = 60.0d-6   ! 60.d-06
    phos(1)   =  2.0d-6   !  2.d-06
    depth(1)  = 0.d0
    Patm(1)   = 1.0d0            !Atmospheric pressure (atm)
    lat(1)    = 0.d0

!  Select input var 'x', choosing set of dy_i/dx to be computed, where y_i are the diff output vars
   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                             &  ! INPUT
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf,                     &
             optGAS='Ppot'    )
   h(1) = 10**(- ph(1))

   write (*,*) "Variables:"
   write (*,*) "          h,           ph,         pco2,         fco2,           co2",&
              "           hco3,           co3,         OmegaA,         OmegaC"
   write (*,"(9ES15.6)")  h, ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC

   write (*,*) "Absolute derivatives" 
   write (*,*) "             dh_dx         dpco2_dx       dfco2_dx         dco2_dx      dhco3_dx       dco3_dx", &
        "       dOmegaA_dx     dOmegaC_dx"

   do i = 1,13
      call derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,                      &
                      dco3_dx, dOmegaA_dx, dOmegaC_dx,                                   &
                      temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1, devar(i),     &
                      optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf          )

      write (*,"(A3,A5,8ES15.6)")  devar(i), "  :  ", dh_dx(1), dpco2_dx(1), dfco2_dx(1), dco2_dx(1), dhco3_dx(1), dco3_dx(1), &
          dOmegaA_dx(1), dOmegaC_dx(1)
   end do


  STOP
END PROGRAM test_derivnum

    
