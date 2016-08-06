!> \file num_deriv.f90
!! \BRIEF 
!> Module with derivnum subroutine - compute numerical derivatives of carbonate system vars 
!! with respect to DIC, Alk, total phosphorus, total silicon, T, S
MODULE mderivnum
CONTAINS
!>    Computes numerically partial derivatives of carbonate system variables 
!!    with respect to one input variable among: DIC, Alk, Silicate, phosphate, temperature and salinity.
!!    Computes derivatives of   (H+, CO2*, HCO3- and CO32-, OmegaA, OmegaC)    as 1D arrays 
!!    Accepts as input the same as subroutine vars(), among which are  : temperature, salinity, pressure,
!!    total alkalinity (ALK), dissolved inorganic carbon (DIC),
!!    silica and phosphate concentrations (all 1-D arrays)
!!    It accepts also a 3-character code identifying the variable with respect to which derivative is requested.
SUBROUTINE derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,                   &
                     dco3_dx, dOmegaA_dx, dOmegaC_dx,                                &
                     temp, sal, alk, dic, sil, phos, Patm, depth, lat, N, derivar,   &
                     optCON, optT, optP, optB, optK1K2, optKf, optGAS                )
                
  ! Computes numerical derivatives by the centered difference method.
  !
  !     INPUT variables:
  !     ================
  !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP='db')
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
  !             = dummy array (unused when optP='db')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not in situ temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !     sal     = salinity in [psu]
  !     alk     = total alkalinity in [eq/m^3] with optCON = 'mol/m3'
  !             =               [eq/kg]  with optCON = 'mol/kg'
  !     dic     = dissolved inorganic carbon [mol/m^3] with optCON = 'mol/m3'
  !             =                            [mol/kg]  with optCON = 'mol/kg'
  !     sil     = silica    [mol/m^3] with optCON = 'mol/m3'
  !             =           [mol/kg]  with optCON = 'mol/kg'
  !     phos    = phosphate [mol/m^3] with optCON = 'mol/m3'
  !             =           [mol/kg]  with optCON = 'mol/kg'
  !     Patm    = atmospheric pressure [atm]
  !     derivar = 3-character identifier of input variable with respect to which derivative is requested
  !               possibilities are 'alk', 'dic', 'pho', 'sil', 'tem', or 'sal'
  !
  !     INPUT options:
  !     ==============
  !     -----------
  !     optCON: choose input & output concentration units - mol/kg (data) vs. mol/m^3 (models)
  !     -----------
  !       -> 'mol/kg' for DIC, ALK, sil, & phos given on mokal scale, i.e., in mol/kg  (std DATA units)
  !       -> 'mol/m3' for DIC, ALK, sil, & phos given in mol/m^3 (std MODEL units)
  !     -----------
  !     optT: choose in situ vs. potential temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
  !       -> 'Tinsitu' means input is already in-situ Temperature, not pot. Temp ("tempis" not computed)
  !     ---------
  !     optP: choose depth (m) vs pressure (db) as input
  !     ---------
  !       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
  !       -> 'db' means "depth" input is already in situ pressure [db], not m (thus p = depth)
  !     ---------
  !     optB: choose total boron formulation - Uppström (1974) vs. Lee et al. (2010)
  !     ---------
  !       -> 'u74' means use classic formulation of Uppström (1974) for total Boron
  !       -> 'l10' means use newer   formulation of Lee et al. (2010) for total Boron
  !     ---------
  !     optK1K2:
  !     ---------
  !       -> 'l'   means use Lueker et al. (2000) formulations for K1 & K2 (recommended by Dickson et al. 2007)
  !                **** BUT this should only be used when 2 < T < 35 and 19 < S < 43
  !       -> 'm10' means use Millero (2010) formulation for K1 & K2 (see Dickson et al., 2007)
  !                **** Valid for 0 < T < 50°C and 1 < S < 50 psu
  !     ----------
  !     optKf:
  !     ----------
  !       -> 'pf' means use Perez & Fraga (1987) formulation for Kf (recommended by Dickson et al., 2007)
  !               **** BUT Valid for  9 < T < 33°C and 10 < S < 40.
  !       -> 'dg' means use Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)
  !     -----------
  !     optGAS: choose in situ vs. potential fCO2 and pCO2
  !     ---------
  !       PRESSURE corrections for K0 and the fugacity coefficient (Cf) 
  !       -> 'Pzero'   = 'zero order' fCO2 and pCO2 (typical approach, which is flawed)
  !                      considers in situ T & only atm pressure (hydrostatic=0)
  !       -> 'Ppot'    = 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !                      considers potential T & only atm pressure (hydrostatic press = 0)
  !       -> 'Pinsitu' = 'in situ' fCO2 and pCO2 (accounts for huge effects of pressure)
  !                      considers in situ T & total pressure (atm + hydrostatic)
  !     ---------

  !     OUTPUT variables:
  !     =================
  !     h    = derivative of H+ concentration on total scale 
  !            h is expressed in [mol/kg] or [mol/m^3] depending on optCON
  !     pco2 = derivative of CO2 partial pressure (uatm/...)
  !     fco2 = derivative of CO2 fugacity (uatm/...)
  !     co2  = derivative of aqueous CO2 concentration in [mol/kg/...] or [mol/m^3/...] depending on optCON
  !     hco3 = derivative of bicarbonate (HCO3-) concentration in [mol/kg/...] or [mol/m^3/...] depending on optCON
  !     co3  = derivative of carbonate (CO3--) concentration in [mol/kg/...] or [mol/m^3/...] depending on optCON
  !     OmegaA = derivative of Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC = derivative of Omega for calcite, i.e., the   calcite saturation state

  USE msingledouble
  USE mvars

  IMPLICIT NONE

! Input variables
  !>     number of records
  INTEGER, INTENT(in) :: N
  !> either <b>in situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) <b>[degree C]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: temp
  !> salinity <b>[psu]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sal
  !> total alkalinity in <b>[eq/m^3]</b> (when optCON = 'mol/m3') OR in <b>[eq/kg]</b>  (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: alk
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: dic
  !> SiO2 concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sil
  !> phosphate concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: phos
!f2py optional , depend(sal) :: n=len(sal)
  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: Patm
  !> depth in \b meters (when optP='m') or \b decibars (when optP='db')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: lat
  ! 3-character identifier of input variable with respect to which derivative is requested
  ! ('alk', 'dic', 'pho', 'sil', 'tem', 'sal'
  CHARACTER(3), INTENT(in) ::  derivar

  !> choose either \b 'mol/kg' (std DATA units) or \b 'mol/m3' (std MODEL units) to select 
  !! concentration units for input (for alk, dic, sil, phos) & output (co2, hco3, co3)
  CHARACTER(6), INTENT(in) :: optCON
  !> choose \b 'Tinsitu' for in situ temperature or \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
  CHARACTER(7), INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
  CHARACTER(2), INTENT(in) :: optP
  !> for total boron, choose either \b 'u74' (Uppstrom, 1974) or \b 'l10' (Lee et al., 2010).
  !! The 'l10' formulation is based on 139 measurements (instead of 20), 
  !! uses a more accurate method, and
  !! generally increases total boron in seawater by 4% 
!f2py character*3 optional, intent(in) :: optB='l10'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optB
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
!f2py character*2 optional, intent(in) :: optKf='pf'
  CHARACTER(2), OPTIONAL, INTENT(in) :: optKf
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) 
!f2py character*3 optional, intent(in) :: optK1K2='l'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optK1K2
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction) 
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optGAS

! Output variables:
  !> derivative of H+ concentration on the <b>total scale</b>
  !! either in <b>[mol/m^3/...]</b> or <b>[mol/kg/...</b>] depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: dh_dx
  !> derivative of CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: dpco2_dx
  !> derivative of CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: dfco2_dx
  !> derivative of aqueous CO2* concentration
  !! either in <b>[mol/m^3/...]</b> or <b>[mol/kg/...</b>] depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: dco2_dx
  !> derivative of bicarbonate ion (HCO3-) concentration 
  !! either in <b>[mol/m^3/..]</b> or <b>[mol/kg/...]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: dhco3_dx
  !> derivative of carbonate ion (CO3--) concentration
  !! either in <b>[mol/m^3/...]</b> or <b>[mol/kg/...]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: dco3_dx
  !> derivative of Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: dOmegaA_dx
  !> derivative of Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: dOmegaC_dx

  ! Local variables

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS

   INTEGER :: i
! Relative delta_x on input variable : constant ratio
  REAL(kind=rx) :: rel_delta_x
! Default value of relative delta_x
#if USE_PRECISION == 2
  REAL(kind=rx), PARAMETER :: def_rel_delta = 1.d-6
#else
  REAL(kind=rx), PARAMETER :: def_rel_delta = 1.e-4
#endif
  
  ! Actual delta_x on input variable
  REAL(kind=rx), DIMENSION(N) :: delta_x
  
  ! Input vars at low and high points
  REAL(kind=rx), DIMENSION(N) :: dic_1, dic_2
  REAL(kind=rx), DIMENSION(N) :: alk_1, alk_2
  REAL(kind=rx), DIMENSION(N) :: sil_1, sil_2
  REAL(kind=rx), DIMENSION(N) :: phos_1, phos_2
  REAL(kind=rx), DIMENSION(N) :: temp_1, temp_2
  REAL(kind=rx), DIMENSION(N) :: sal_1, sal_2
  
  ! Output carbonate vars at low and high points
  REAL(kind=rx), DIMENSION(N) :: ph_1, ph_2, H_1, H_2
  REAL(kind=rx), DIMENSION(N) :: pco2_1, pco2_2
  REAL(kind=rx), DIMENSION(N) :: fco2_1, fco2_2
  REAL(kind=rx), DIMENSION(N) :: co2_1, co2_2
  REAL(kind=rx), DIMENSION(N) :: hco3_1, hco3_2
  REAL(kind=rx), DIMENSION(N) :: co3_1, co3_2
  REAL(kind=rx), DIMENSION(N) :: OmegaA_1, OmegaA_2
  REAL(kind=rx), DIMENSION(N) :: OmegaC_1, OmegaC_2

  !> Revelle factor, i.e., dpCO2/pCO2 / dDIC/DIC
  REAL(kind=rx), DIMENSION(N) :: BetaD
  !> in-situ density of seawater; rhoSW = f(s, t, p) in <b>[kg/m3]</b>
  REAL(kind=rx), DIMENSION(N) :: rhoSW
  !> pressure <b>[decibars]</b>; p = f(depth, latitude) if computed from depth [m] (when optP='m') OR p = depth [db] (when optP='db')
  REAL(kind=rx), DIMENSION(N) :: p
  !> in-situ temperature \b <b>[degrees C]</b>
  REAL(kind=rx), DIMENSION(N) :: tempis
  
!----------------------------------------------------------------------------------

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaraion
  IF (PRESENT(optB)) THEN
!   print *,"optB present:"
!   print *,"optB = ", optB 
    opB = optB
  ELSE
!   Default is Lee et al (2010) for total boron
!   print *,"optB NOT present:"
    opB = 'l10'
!   print *,"opB = ", opB 
  ENDIF
  IF (PRESENT(optKf)) THEN
!   print *,"optKf = ", optKf
    opKf = optKf
  ELSE
!   print *,"optKf NOT present:"
!   Default is Perez & Fraga (1987) for Kf
    opKf = 'pf'
!   print *,"opKf = ", opKf
  ENDIF
  IF (PRESENT(optK1K2)) THEN
!   print *,"optK1K2 = ", optK1K2
    opK1K2 = optK1K2
  ELSE
!   print *,"optK1K2 NOT present:"
!   Default is Lueker et al. 2000) for K1 & K2
    opK1K2 = 'l'
!   print *,"opK1K2 = ", opK1K2
  ENDIF
  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF

  ! Determine HICH and LOW points
  ! -----------------------------

  ! Default values : not perturbed
  ! no change on DIC and Alk
  dic_1(:) = dic(:)
  dic_2(:) = dic(:)
  alk_1(:) = alk(:)
  alk_2(:) = alk(:)
  ! no change on Sil total and Phos total
  sil_1(:) = sil(:)
  sil_2(:) = sil(:)
  phos_1(:) = phos(:)
  phos_2(:) = phos(:)
  ! no change on T and S
  temp_1(:) = temp(:)
  temp_2(:) = temp(:)
  sal_1(:) = sal(:)
  sal_2(:) = sal(:)
  
  ! Choose relative delta_x
  SELECT CASE (derivar)
      CASE ('alk')
          rel_delta_x = def_rel_delta
      CASE ('dic')
          rel_delta_x = def_rel_delta
      CASE ('pho')
          rel_delta_x = 1.0d-3
      CASE ('sil')
          rel_delta_x = 1.0d-3
      CASE ('tem')
          rel_delta_x = 1.0d-4
      CASE ('sal')
          rel_delta_x = 1.0d-4
      CASE DEFAULT
          PRINT *,"derivar must be 3-char variable: 'alk', 'dic', 'pho', 'sil', 'tem', or 'sal'"
          STOP
  END SELECT

  ! Compute HIGH and LOW points
  SELECT CASE (derivar)
    CASE ('alk')
      ! change slightly Alk
      DO i = 1, N
        alk_1(i) = alk(i) - alk(i) * rel_delta_x
        alk_2(i) = alk(i) + alk(i) * rel_delta_x
        delta_x(i) = alk_2(i) - alk_1(i)
      END DO

    CASE ('dic')
      ! change slightly DIC
      DO i = 1, N
        dic_1(i) = dic(i) - dic(i) * rel_delta_x
        dic_2(i) = dic(i) + dic(i) * rel_delta_x
        delta_x(i) = dic_2(i) - dic_1(i)
      END DO

    CASE ('sil')
      ! change slightly Silicate concentration
      DO i = 1, N
        sil_1(i) = sil(i) - sil(i) * rel_delta_x
        sil_2(i) = sil(i) + sil(i) * rel_delta_x
        delta_x(i) = sil_2(i) - sil_1(i)
      END DO

    CASE ('pho')
      ! change slightly Phosphate concentration
      DO i = 1, N
        phos_1(i) = phos(i) - phos(i) * rel_delta_x
        phos_2(i) = phos(i) + phos(i) * rel_delta_x
        delta_x(i) = phos_2(i) - phos_1(i)
      END DO

    CASE ('tem')
      ! change slightly Temperature
      DO i = 1, N
        temp_1(i) = temp(i) - temp(i) * rel_delta_x
        temp_2(i) = temp(i) + temp(i) * rel_delta_x
        delta_x(i) = temp_2(i) - temp_1(i)
      END DO

    CASE ('sal')
      ! change slightly the salinity
      DO i = 1, N
        sal_1(i) = sal(i) - sal(i) * rel_delta_x
        sal_2(i) = sal(i) + sal(i) * rel_delta_x
        delta_x(i) = sal_2(i) - sal_1(i)
      END DO

  END SELECT
  
  ! PERTURBATION:
  !-------------
  ! Point 1: (one of Alk, DIC, sil, phos, T or S is somewhat smaller)
  call vars(ph_1, pco2_1, fco2_1, co2_1, hco3_1, co3_1, OmegaA_1, OmegaC_1, BetaD, rhoSW, p, tempis,  &  ! OUTPUT
        temp_1, sal_1, alk_1, dic_1, sil_1, phos_1, Patm, depth, lat, N,                      &  ! INPUT
        optCON, optT, optP, optB, optK1K2, optKf)                                    ! INPUT OPTIONS

  ! Point 2: (one of Alk, DIC, sil, phos, T or S is somewhat bigger)
  call vars(ph_2, pco2_2, fco2_2, co2_2, hco3_2, co3_2, OmegaA_2, OmegaC_2, BetaD, rhoSW, p, tempis,  &  ! OUTPUT
        temp_2, sal_2, alk_2, dic_2, sil_2, phos_2, Patm, depth, lat, N,                      &  ! INPUT
        optCON, optT, optP, optB, optK1K2, optKf)                                    ! INPUT OPTIONS

  ! Compute [H+] in mol/Kg
  DO i = 1, N
    H_1(i) = 10.**(-ph_1(i))
    H_2(i) = 10.**(-ph_2(i))
  END DO

  DO i = 1, N
    ! Compute ratio dy/dx with dy = Centered difference
    dh_dx(i)    = (H_2(i)    - H_1(i))    / delta_x(i)
    dpco2_dx(i) = (pco2_2(i) - pco2_1(i)) / delta_x(i)
    dfco2_dx(i) = (fco2_2(i) - fco2_1(i)) / delta_x(i)
    dco2_dx(i)  = (co2_2(i)  - co2_1(i))  / delta_x(i)
    dhco3_dx(i) = (hco3_2(i) - hco3_1(i)) / delta_x(i)
    dco3_dx(i)  = (co3_2(i)  - co3_1(i))  / delta_x(i)
    dOmegaA_dx(i) = (OmegaA_2(i) - OmegaA_1(i)) / delta_x(i)
    dOmegaC_dx(i) = (OmegaC_2(i) - OmegaC_1(i)) / delta_x(i)
  END DO
  
  IF  (trim(optCON) == 'mol/m3') THEN
    ! convert dH/dx from mol/Kg/... to mol/m3/...
    DO i = 1, N
      dh_dx(i) = dh_dx(i) * rhoSW(i)
    END DO
  ENDIF
  
  
  RETURN
END SUBROUTINE derivnum
END MODULE mderivnum
