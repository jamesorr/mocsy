!> \file vars.f90
!! \BRIEF 
!> Module with vars subroutine - compute carbonate system vars from DIC,Alk,T,S,P,nuts
MODULE mvars
CONTAINS
!>    Computes standard carbonate system variables (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
!!    as 1D arrays FROM
!!    temperature, salinity, pressure,
!!    total alkalinity (ALK), dissolved inorganic carbon (DIC),
!!    silica and phosphate concentrations (all 1-D arrays)
SUBROUTINE vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                      &
                optCON, optT, optP, optB, optK1K2, optKf, optGAS                          )

  !   Purpose:
  !     Computes other standard carbonate system variables (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
  !     as 1D arrays
  !     FROM:
  !     temperature, salinity, pressure,
  !     total alkalinity (ALK), dissolved inorganic carbon (DIC),
  !     silica and phosphate concentrations (all 1-D arrays)

  !     INPUT variables:
  !     ================
  !     Patm    = atmospheric pressure [atm]
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
  !     ph   = pH on total scale
  !     pco2 = CO2 partial pressure (uatm)
  !     fco2 = CO2 fugacity (uatm)
  !     co2  = aqueous CO2 concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     hco3 = bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     co3  = carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC = Omega for calcite, i.e., the   calcite saturation state
  !     BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
  !     rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
  !     p = pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
  !     tempis  = in-situ temperature [degrees C]

  USE msingledouble
  USE mconstants
  USE mp80
  USE mrho
  USE msw_temp
  USE mvarsolver

  IMPLICIT NONE

! Input variables
  !>     number of records
  INTEGER, INTENT(in) :: N
  !> either <b>in situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) <b>[degree C]</b>
  REAL(kind=r4), INTENT(in),    DIMENSION(N) :: temp
  !> salinity <b>[psu]</b>
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: sal
  !> total alkalinity in <b>[eq/m^3]</b> (when optCON = 'mol/m3') OR in <b>[eq/kg]</b>  (when optCON = 'mol/kg')
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: alk
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: dic
  !> SiO2 concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: sil
  !> phosphate concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: phos
!f2py optional , depend(sal) :: n=len(sal)
  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: Patm
  !> depth in \b meters (when optP='m') or \b decibars (when optP='db')
  REAL(kind=r4), INTENT(in),    DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=r4), INTENT(in),    DIMENSION(N) :: lat

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
  !> pH on the <b>total scale</b>
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: ph
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: OmegaC
  !> Revelle factor, i.e., dpCO2/pCO2 / dDIC/DIC
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: BetaD
  !> in-situ density of seawater; rhoSW = f(s, t, p) in <b>[kg/m3]</b>
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: rhoSW
  !> pressure <b>[decibars]</b>; p = f(depth, latitude) if computed from depth [m] (when optP='m') OR p = depth [db] (when optP='db')
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: p
  !> in-situ temperature \b <b>[degrees C]</b>
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: tempis

! Local variables
  REAL(kind=r4) :: ssal, salk, sdic, ssil, sphos

  REAL(kind=r4) :: tempot, tempis68, tempot68
! REAL(kind=r8) :: dtempot, dtempot68
  REAL(kind=r8) :: drho

  REAL(kind=r8) :: K0, K1, K2, Kb, Kw, Ks, Kf, Kspc
  REAL(kind=r8) :: Kspa, K1p, K2p, K3p, Ksi
  REAL(kind=r8) :: St, Ft, Bt

  REAL(kind=r8), DIMENSION(1) :: aK0, aK1, aK2, aKb, aKw, aKs, aKf, aKspc
  REAL(kind=r8), DIMENSION(1) :: aKspa, aK1p, aK2p, aK3p, aKsi
  REAL(kind=r8), DIMENSION(1) :: aSt, aFt, aBt

  REAL(kind=r8) :: Patmd, Ptot, Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  REAL(kind=r8) :: Phydro_atm

  INTEGER :: i, icount

  REAL(kind=r8) :: t, tk, prb
  REAL(kind=r8) :: s
  REAL(kind=r8) :: tc, ta
  REAL(kind=r8) :: sit, pt
  REAL(kind=r8) :: Hinit
  REAL(kind=r8) :: ah1

  REAL(kind=r8) :: HSO4, HF, HSI, HPO4
  REAL(kind=r8) :: ab, aw, ac, ah2, erel

  REAL(kind=r8) :: cu, cb, cc

  REAL(kind=r8), DIMENSION(2) :: dicdel, pco2del
  REAL(kind=r8) :: dx, Rf
  REAL(kind=r8) :: dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC

  INTEGER :: kcomp
  INTEGER :: j, minusplus

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS

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

  icount = 0
  DO i = 1, N
     icount = icount + 1
!    ===============================================================
!    Convert model depth -> press; convert model Theta -> T in situ
!    ===============================================================
!    * Model temperature tracer is usually "potential temperature"
!    * Model vertical grid is usually in meters
!    BUT carbonate chem routines require pressure & in-situ T
!    Thus before computing chemistry, if appropriate,
!    convert these 2 model vars (input to this routine)
!    - depth [m] => convert to pressure [db]
!    - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
     !print *,"optP =", optP, "end"
     IF (trim(optP) == 'm' ) THEN
!       Compute pressure [db] from depth [m] and latitude [degrees]
        p(i) = p80(depth(i), lat(i))
     ELSEIF (trim(optP) == 'db') THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        p(i) = depth(i)
     ELSE
        !print *,"optP =", optP, "end"
        PRINT *,"optP must be either 'm' or 'db'"
        STOP
     ENDIF

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
     IF (trim(optT) == 'Tpot' .OR. trim(optT) == 'tpot') THEN
        tempot = temp(i)
!       This is the case for most models and some data
!       a) Convert the pot. temp on today's "ITS 90" scale to older IPTS 68 scale
!          (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
        tempot68 = (tempot - 0.0002) / 0.99975
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
        tempis68 = sw_temp(sal(i), tempot68, p(i), 0. )
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis(i) = 0.99975*tempis68 + 0.0002
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (trim(optT) == 'Tinsitu' .OR. trim(optT) == 'tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis(i) = temp(i)
        tempis68  = (temp(i) - 0.0002) / 0.99975
!       dtempot68 = sw_ptmp(DBLE(sal(i)), DBLE(tempis68), DBLE(p), 0.0d0)
!       dtempot   = 0.99975*dtempot68 + 0.0002
     ELSE
        PRINT *,"optT must be either 'Tpot' or 'Tinsitu'"
        PRINT *,"you specified optT =", trim(optT) 
        STOP
     ENDIF

!    ================================================================
!    Carbonate chemistry computations
!    ================================================================
     IF (dic(i) > 0. .AND. dic(i) < 1.0e+4) THEN
!       Test to indicate if any of input variables are unreasonable
        IF (       sal(i) < 0.   &
             .OR.  alk(i) < 0.   &
             .OR.  dic(i) < 0.   &
             .OR.  sil(i) < 0.   &
             .OR. phos(i) < 0.   &
             .OR.  sal(i) > 1e+3 &
             .OR.  alk(i) > 1e+3 &
             .OR.  dic(i) > 1e+3 &
             .OR.  sil(i) > 1e+3 &
             .OR. phos(i) > 1e+3) THEN
           PRINT *, 'i, icount, tempot, sal,    alk,    dic,    sil,    phos =', &
                     i, icount, tempot, sal(i), alk(i), dic(i), sil(i), phos(i)
        ENDIF
!       Zero out any negative salinity, phosphate, silica, dic, and alk
        IF (sal(i) < 0.0) THEN
           ssal = 0.0
        ELSE
           ssal = sal(i)
        ENDIF
        IF (phos(i) < 0.0) THEN
           sphos = 0.0
        ELSE
           sphos = phos(i)
        ENDIF
        IF (sil(i) < 0.0) THEN
           ssil = 0.0
        ELSE
           ssil = sil(i)
        ENDIF
        IF (dic(i) < 0.0) THEN
          sdic = 0.0
        ELSE
          sdic = dic(i)
        ENDIF
        IF (alk(i) < 0.0) THEN
          salk = 0.0
        ELSE
          salk = alk(i)
        ENDIF

!       Absolute temperature (Kelvin) & related variables
        t  = DBLE(tempis(i))
        tk = 273.15d0 + t

!       Atmospheric pressure
        Patmd = DBLE(Patm(i))
!       Hydrostatic pressure (prb is in bars)
        prb = DBLE(p(i)) / 10.0d0
        Phydro_atm = prb / 1.01325d0  ! convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
!       Total pressure [atm]
        IF     (trim(opGAS) == 'Pzero'     .OR. trim(opGAS) == 'pzero') THEN
           Ptot = Patmd               ! total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Ppot'    .OR. trim(opGAS) == 'ppot') THEN
           Ptot = Patmd               ! total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
           Ptot = Patmd + Phydro_atm   ! total pressure (in atm) = atmospheric pressure + hydrostatic pressure
        ELSE
           PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
           STOP
        ENDIF

!       Salinity (equivalent array in double precision)
        s = DBLE(ssal)

!       Get all equilibrium constants and total concentrations of SO4, F, B
        CALL constants(aK0, aK1, aK2, aKb, aKw, aKs, aKf, aKspc, aKspa,  &
                       aK1p, aK2p, aK3p, aKsi,                           &
                       aSt, aFt, aBt,                                    &
                       temp(i), sal(i), Patm(i),                         &
                       depth(i), lat(i), 1,                              &
                       optT, optP, opB, opK1K2, opKf, opGAS              )

!       Unlike f77, in F90 we can't assign an array (dimen=1) to a scalar in a routine argument
!       Thus, set scalar constants equal to array (dimension=1) values required as arguments
        K0 = aK0(1) ; K1 = aK1(1) ; K2 = aK2(1) ; Kb = aKb(1) ; Kw = aKw(1) 
        Ks = aKs(1) ; Kf = aKf(1) ; Kspc = aKspc(1) ; Kspa = aKspa(1) 
        K1p = aK1p(1) ; K2p = aK2p(1) ; K3p = aK3p(1) ; Ksi = aKsi(1)
        St = aSt(1) ; Ft = aFt(1) ; Bt = aBt(1)

!       Compute in-situ density [kg/m^3]
        rhoSW(i) =  rho(ssal, tempis68, REAL(prb))

!       Either convert units of DIC and ALK (MODEL case) or not (DATA case)
        IF     (trim(optCON) == 'mol/kg') THEN
!          No conversion:
!          print *,'DIC and ALK already given in mol/kg (std DATA units)'
           drho = 1.
        ELSEIF (trim(optCON) == 'mol/m3') THEN
!          Do conversion:
!          print *,"DIC and ALK given in mol/m^3 (std MODEL units)"
           drho = DBLE(rhoSW(i))
        ELSE
           PRINT *,"optCON must be either 'mol/kg' or 'mol/m3'"
           STOP
        ENDIF

        tc  = DBLE(sdic)/drho
        ta  = DBLE(salk)/drho
        sit = DBLE(ssil)/drho
        pt  = DBLE(sphos)/drho

!       Solve for pH and all other variables
!       ------------------------------------
        CALL varsolver(dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC,     &
                      t, s, ta, tc, pt, sit,                                       &
                      Bt, St, Ft,                                                  &
                      K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi,  & 
                      Patmd, prb, drho, opGAS                                     )

!       Convert all output variables from double to single precision
        pH(i)     = REAL(dph)
        co2(i)    = REAL(dco2)
        hco3(i)   = REAL(dhco3)
        co3(i)    = REAL(dco3)
        fCO2(i)   = REAL(dfCO2)
        pCO2(i)   = REAL(dpCO2)
        OmegaA(i) = REAL(dOmegaA)
        OmegaC(i) = REAL(dOmegaC)

!       Compute Revelle factor numerically (derivative using centered-difference scheme)
        DO j=1,2
           minusplus = (-1)**j
           dx = 0.1 * 1e-6         ! Numerical tests found for DIC that optimal dx = 0.1 umol/kg (0.1e-6 mol/kg)
           dicdel(j) = tc + DBLE(minusplus)*dx/2.0d0
            CALL varsolver(dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC, &
               t, s, ta, dicdel(j), pt, sit,                                       &
               Bt, St, Ft,                                                         &
               K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi,         & 
               Patmd, prb, drho, optGAS                                            )
            pco2del(j) = dpco2
        END DO
       !Classic finite centered difference formula for derivative (2nd order accurate)
        Rf = (pco2del(2) - pco2del(1)) / (dicdel(2) - dicdel(1))       ! dpCO2/dDIC
       !Rf = (pco2del(2) - pco2del(1)) / (dx)                          ! dpCO2/dDIC (same as just above)
        Rf = Rf * tc / dpco2                                           ! R = (dpCO2/dDIC) * (DIC/pCO2)

        BetaD(i) = REAL(Rf)

     ELSE

        ph(i)     = 1.e20_r4
        pco2(i)   = 1.e20_r4
        fco2(i)   = 1.e20_r4
        co2(i)    = 1.e20_r4
        hco3(i)   = 1.e20_r4
        co3(i)    = 1.e20_r4
        OmegaA(i) = 1.e20_r4
        OmegaC(i) = 1.e20_r4
        BetaD(i)  = 1.e20_r4
        rhoSW(i)  = 1.e20_r4
        p(i)      = 1.e20_r4
        tempis(i) = 1.e20_r4

     ENDIF

  END DO

  RETURN
END SUBROUTINE vars
END MODULE mvars
