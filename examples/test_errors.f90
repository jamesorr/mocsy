!> \file test_errors.f90
!! \BRIEF Fortran 90 program to test errors.f90
PROGRAM test_errors

   USE msingledouble
   USE merrors

   IMPLICIT NONE


!  Output variables:
   REAL(kind=rx), DIMENSION(6) :: eh, epco2, efco2, eco2, ehco3, eco3, eOmegaA, eOmegaC
!  Input variables
   REAL(kind=rx), DIMENSION(6) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
   REAL(kind=rx), DIMENSION(6) :: temp_e, sal_e, ALK_e, DIC_e, sil_e, phos_e
   REAL(kind=rx), DIMENSION(7) :: epK, epKstd, epK0
   REAL(kind=rx)               :: ebt
   
!  Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2

!  Local variables:
   INTEGER ::  i

!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'm'       ! input variable 'depth' is in meters
   optB    = 'u74'
   optK1K2 = 'l'
   optKf   = 'dg'

    ! Input errors
    ALK_e(:)  = 2.d-6    ! (2 umol/kg for ALK)
    DIC_e(:)  = 2.d-6    ! (2 umol/kg for DIC)
    sal_e(:)  = 0.0d0     ! (psu)
    temp_e(:) = 0.0d0    ! (C)
    phos_e(:) = 0.1d-6
    sil_e(:)  = 4.0d-6    
 !  epK(:)    = 0.0
 !  epK = (/0.004_r8, 0.015_r8, 0.03_r8, 0.01_r8, 0.01_r8, 0.02_r8, 0.02_r8/)
    epKstd = (/0.004_r8, 0.015_r8, 0.03_r8, 0.01_r8, 0.01_r8, 0.02_r8, 0.02_r8/)
    epK0 = epKstd * 0.0d0
    ebt = 0.01_rx

   ! ------------------
   ! 1s test : 1 record 
   ! ------------------
    
    temp(1)   = 18.0d0           ! Can be "Potential temperature" or "In situ temperature" (see optT below)
    sal(1)    = 35.0d0           ! Salinity (practical scale)
    alk(1)    = 2300.d-6         ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
    dic(1)    = 2000.d-6         ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
!
    sil(1)    = 60.d-6           ! 60
    phos(1)   = 2.d-6            !  2
!
!   sil(1)    = 0.d0
!   phos(1)   = 0.d0
!
    depth(1)  = 0.d0
    Patm(1)   = 1.0d0            ! Atmospheric pressure (atm)
    lat(1)    = 0.d0

!  Compute output errors

!  -----------------------------------------------------------------------------------------------------------------------------
   write (*,*) "Test 1: NO error for constants nor for Bt"
!  -----------------------------------------------------------------------------------------------------------------------------
   call errors(eh, epco2, efco2, eco2, ehco3, eco3, eOmegaA, eOmegaC,            &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                &  ! INPUT
             temp_e, sal_e, ALK_e, DIC_e, sil_e, phos_e,                         &
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf, epK=epK0, ebt=0.0_rx )
   write (*,*) "              eh          epco2          efco2          eco2           ehco3          eco3   ", &
        "      eOmegaA        eOmegaC   "
   DO i = 1,1
     write (*,"(5x,8ES15.6)")  eh(i), epco2(i), efco2(i), eco2(i), ehco3(i), eco3(i), eOmegaA(i), eOmegaC(i)
   END DO
   write (*,*) " "

!  -----------------------------------------------------------------------------------------------------------------------------
   write (*,*) "Test 2: Default error on constants but Bt=0"
!  -----------------------------------------------------------------------------------------------------------------------------
   call errors(eh, epco2, efco2, eco2, ehco3, eco3, eOmegaA, eOmegaC,          &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,              &  ! INPUT
             temp_e, sal_e, ALK_e, DIC_e, sil_e, phos_e,                       &
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf, epK=epKstd, ebt=0.0_rx )
   write (*,*) "              eh          epco2          efco2          eco2           ehco3          eco3   ", &
        "      eOmegaA        eOmegaC   "
   DO i = 1,1
     write (*,"(5x,8ES15.6)")  eh(i), epco2(i), efco2(i), eco2(i), ehco3(i), eco3(i), eOmegaA(i), eOmegaC(i)
   END DO
   write (*,*) " "

!  -----------------------------------------------------------------------------------------------------------------------------
   write (*,*) "Test 3: Default error on constants and Bt (specified as arguments) - ebt = 0.01"
!  -----------------------------------------------------------------------------------------------------------------------------
   call errors(eh, epco2, efco2, eco2, ehco3, eco3, eOmegaA, eOmegaC,          &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,              &  ! INPUT
             temp_e, sal_e, ALK_e, DIC_e, sil_e, phos_e,                       &
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf, epK=epKstd, ebt=ebt )
   write (*,*) "              eh          epco2          efco2          eco2           ehco3          eco3   ", &
        "      eOmegaA        eOmegaC   "
   DO i = 1,1
     write (*,"(5x,8ES15.6)")  eh(i), epco2(i), efco2(i), eco2(i), ehco3(i), eco3(i), eOmegaA(i), eOmegaC(i)
   END DO
   write (*,*) " "

!  -----------------------------------------------------------------------------------------------------------------------------
   write (*,*) "Test 4: Default error on constants but ebt=0.04"
!  -----------------------------------------------------------------------------------------------------------------------------
   call errors(eh, epco2, efco2, eco2, ehco3, eco3, eOmegaA, eOmegaC,          &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,              &  ! INPUT
             temp_e, sal_e, ALK_e, DIC_e, sil_e, phos_e,                       &
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf, epK=epKstd, ebt=0.04_rx )
   write (*,*) "              eh          epco2          efco2          eco2           ehco3          eco3   ", &
        "      eOmegaA        eOmegaC   "
   DO i = 1,1
     write (*,"(5x,8ES15.6)")  eh(i), epco2(i), efco2(i), eco2(i), ehco3(i), eco3(i), eOmegaA(i), eOmegaC(i)
   END DO
   write (*,*) " "

!  -----------------------------------------------------------------------------------------------------------------------------
   write (*,*) "Test 5: Default error on constants and Bt (not specified as arguments) - should give same results aas Test 3"
!  -----------------------------------------------------------------------------------------------------------------------------
   call errors(eh, epco2, efco2, eco2, ehco3, eco3, eOmegaA, eOmegaC,          &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,              &  ! INPUT
             temp_e, sal_e, ALK_e, DIC_e, sil_e, phos_e,                       &
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf       )
   write (*,*) "              eh          epco2          efco2          eco2           ehco3          eco3   ", &
        "      eOmegaA        eOmegaC   "
   DO i = 1,1
     write (*,"(5x,8ES15.6)")  eh(i), epco2(i), efco2(i), eco2(i), ehco3(i), eco3(i), eOmegaA(i), eOmegaC(i)
   END DO
   write (*,*) " "

   ! ----------------------------------------
   ! Test 6: six records at increasing depth
   ! ----------------------------------------
    
   ! Simple input data (with CONCENTRATION units typical for DATA)
   ! (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
   DO i = 1,6
     temp(i)   = 2.0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
     sal(i)    = 35.0           !Salinity (practical scale)
     alk(i)    = 2295.*1.e-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
     dic(i)    = 2154.*1.e-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
     sil(i)    = 0.
     phos(i)   = 0.
     depth(i) = real(i-1) * 1000. ! Vary depth from 0 to 5000 db by 1000 db
     Patm(i)   = 1.0            !Atmospheric pressure (atm)
     lat(i)    = 0.
   END DO
   
   ! compute output errors
   call errors(eh, epco2, efco2, eco2, ehco3, eco3, eOmegaA, eOmegaC,          &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 6,              &  ! INPUT
             temp_e, sal_e, ALK_e, DIC_e, sil_e, phos_e,                       &
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf   )

   write (*,*) "Test 6: six records - same input (Southern Ocean surface averages) except depth increases from 0 to 5000 m"
   write (*,*) "depth         eh          epco2          efco2          eco2           ehco3          eco3   ", &
        "      eOmegaA        eOmegaC   "

   DO i = 1,6
     write (*,"(f5.0, 8ES15.6)")  depth(i), eh(i), epco2(i), efco2(i), eco2(i), ehco3(i), eco3(i), eOmegaA(i), eOmegaC(i)
   END DO


  STOP
END PROGRAM test_errors
    
