!> \file gasx.f90
!! \BRIEF 
!> Module with routines needed to compute gas exchange (flxco2, scco2, atmospheric xCO2 and pCO2)
MODULE gasx
CONTAINS
!>    Computes air-sea CO2 flux & surface-ocean carbonate system vars (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
!!    from T, S, P, ALK, DIC, total inorganic silicon, total inorganic phosphorus, all as 1-D arrays
SUBROUTINE flxco2(co2flux, co2ex, dpco2,                                                    &
                  ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
                  temp, sal, alk, dic, sil, phos, kw660, xco2, Patm, dz1, N,                &
                  optCON, optT, optP, optB, optK1K2, optKf, optGAS                          )

  !   Purpose:
  !     Computes air-sea CO2 flux & surface ocean carbonate system vars (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
  !     from T, S, P, ALK, DIC, total inorganic silicon, total inorganic phosphorus, all as 1-D arrays
  !
  !     INPUT variables:
  !     ================
  !     kw660   = gas transfer velocity (piston velocity) for CO2 [m/s] 
  !               without T-dependant Schmidt number correction
  !               but accounting for sea ice fraction (See OCMIP2 design & HOWTO documents for details)
  !     xco2    = atmospheric mole fraction of CO2 [ppm]
  !     Patm    = atmospheric pressure at surface [atm]
  !     dz1     = depth of the top vertical layer of the model [m]
  !     temp    = surface potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not in situ temp)
  !             = surface in situ temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !     sal     = surface salinity in [psu]
  !     alk     = surface total alkalinity in [eq/m^3] with optCON = 'mol/m3'
  !             =                             [eq/kg]  with optCON = 'mol/kg'
  !     dic     = surface dissolved inorganic carbon [mol/m^3] with optCON = 'mol/m3'
  !             =                            [mol/kg]  with optCON = 'mol/kg'
  !     sil     = surface total inorganic silicon   [mol/m^3] with optCON = 'mol/m3'
  !             =                                   [mol/kg]  with optCON = 'mol/kg'
  !     phos    = surface total inorganinc phosphorus [mol/m^3] with optCON = 'mol/m3'
  !             =                                     [mol/kg]  with optCON = 'mol/kg'
  !     INPUT options:
  !     ==============
  !     -----------
  !     optCON: choose input concentration units - mol/kg (data) vs. mol/m^3 (models)
  !     -----------
  !       -> 'mol/kg' for DIC and ALK given on mokal scale, i.e., in mol/kg  (std DATA units)
  !       -> 'mol/m3' for DIC and ALK given in mol/m^3 (std MODEL units)
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
  !                **** Orr et al. (GMDD, 2014) identify large discrepancies between packages w/ this option
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
  !     co2flux = air-to-sea flux of CO2 [mol/(m^2 * s)]
  !     co2ex = time rate of change of surface CO2 due to gas exchange [mol/(m^3 * s)]
  !     dpco2 = difference of oceanic pCO2 minus atmospheric pCO2 [uatm]
  !     ph   = pH on total scale
  !     pco2 = oceanic partial pressure of CO2 (uatm)
  !     fco2 = oceanic fugacity of CO2 (uatm)
  !     co2  = aqueous CO2 concentration [mol/m^3]
  !     hco3 = bicarbonate (HCO3-) concentration [mol/m^3]
  !     co3  = carbonate (CO3--) concentration [mol/m^3]
  !     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC = Omega for calcite, i.e., the   calcite saturation state
  !     BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
  !     rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
  !     p = pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
  !     tempis  = in-situ temperature [degrees C]

  USE msingledouble
  USE mvars
  USE mp2fCO2

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
  !> gas transfer velocity (piston velocity) at a Schmidt number of 660 <b>[m/s]</b>
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: kw660
  !> atmospheric mole fraction of CO2 <b>[ppm]</b>
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: xco2
  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: Patm
  !> thickness of the surface layer of the model <b>[m]</b> 
  REAL(kind=r4), INTENT(in) :: dz1

  !> choose either \b 'mol/kg' (std DATA units) or \b 'mol/m3' (std MODEL units) to select 
  !! concentration units for input (for alk, dic, sil, phos) & output (co2, hco3, co3)
!f2py character*6 optional, intent(in) :: optCON='mol/m3'
  CHARACTER(6), OPTIONAL, INTENT(in) :: optCON
  !> choose \b 'Tinsitu' for in situ temperature or \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
!f2py character*7 optional, intent(in) :: optT='Tinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
!f2py character*2 optional, intent(in) :: optP='m'
  CHARACTER(2), OPTIONAL, INTENT(in) :: optP
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
  !> air-to-sea CO2 flux <b>[mol/(m^2 * s)]</b>
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: co2flux
  !> rate of change of surface DIC concentration <b>[mol/(m^3 * s)]</b>
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: co2ex
  !> difference of surface ocean pCO2 minus atmospheric pCO2 <b>[uatm]</b>
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: dpCO2
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
  REAL(kind=r8) :: tk, invtk, dtemp
  REAL(kind=r8) :: tmp, K0, co2star, co2starair, kwco2
  REAL(kind=r4), DIMENSION(N) :: pCO2atm, fCO2atm
  REAL(kind=r4), DIMENSION(N) :: depth0, lat0
 
  INTEGER :: i
  INTEGER :: kcomp

! Optional arguments: if not present use defaults (defined below)
  CHARACTER(6) :: opCON
  CHARACTER(7) :: opT
  CHARACTER(2) :: opP
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaration
!        Those !f2py statements should be consistent w/ defaults below
  IF (PRESENT(optCON)) THEN
    opCON = optCON
  ELSE
!   Default (typical for models, not data)
    opCON = 'mol/m3'
  ENDIF
  IF (PRESENT(optT)) THEN
    opT = optT
  ELSE
!   Default (this option is irrelevant for surface values, 0 m)
    opT = 'Tinsitu'
  ENDIF
  IF (PRESENT(optP)) THEN
    opP = optP
  ELSE
!   Default (this option is irrelevant for surface values, 0 m)
    opP = 'm'
  ENDIF
  IF (PRESENT(optB)) THEN
    opB = optB
  ELSE
!   Default is Lee et al (2010) for total boron
    opB = 'l10'
  ENDIF
  IF (PRESENT(optKf)) THEN
    opKf = optKf
  ELSE
!   Default is Perez & Fraga (1987) for Kf
    opKf = 'pf'
  ENDIF
  IF (PRESENT(optK1K2)) THEN
    opK1K2 = optK1K2
  ELSE
!   Default is Lueker et al. 2000) for K1 & K2
    opK1K2 = 'l'
  ENDIF
  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF

  depth0 = fco2 * 0.0
  lat0   = depth0
! Compute derived variables from input (DIC, ALK, ...)
  CALL vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
            temp, sal, alk, dic, sil, phos, Patm, depth0, lat0, N,                    &
            opCON, opT, opP, opB, opK1K2, opKf, opGAS                                 )

! Compute pCO2atm [uatm] from xCO2 [ppm], atmospheric pressure [atm], & vapor pressure of seawater
! pCO2atm = (Patm - pH20(i)) * xCO2,   where pH20 is the vapor pressure of seawater [atm]
  CALL x2pCO2atm(xco2, temp, sal, Patm, N, pco2atm)

! Compute fCO2atm [uatm] from pCO2atm [uatm] & fugacity coefficient [unitless]
! fCO2atm = pCO2atm * fugcoeff,   where fugcoeff= exp(Patm*(B + 2.0*xc2*Del)/(R*tk) )
  CALL p2fCO2(pco2atm, temp, Patm, depth0, N, fco2atm)

! Compute flux, absolute rate of change of surface DIC, & Delta pCO2
  DO i = 1, N
     IF (co3(i) .eq. 1.e20_r4) THEN
!       Masked values (land)
        co2flux(i) = 1.e20_r4
        co2ex(i)   = 1.e20_r4
        dpco2(i)   = 1.e20_r4
     ELSE
        dtemp = DBLE(temp(i))
        tk = dtemp + 273.15d0
        invtk = 1.0d0/tk

!       Transfer velocity for CO2 in m/s (see equation [4] in OCMIP2 design document & OCMIP2 Abiotic HOWTO)
        kwco2 = DBLE(kw660(i)) * (660/scco2(dtemp))**0.5

!       Surface K0 [(mol/kg) / atm] at T, S of surface water
        tmp = 9345.17d0*invtk - 60.2409d0 + 23.3585d0 * LOG(tk/100.0d0)
        K0 = EXP( tmp + DBLE(sal(i))*(0.023517d0 - 0.00023656d0*tk + 0.0047036e-4_r8*tk*tk) )

!       "Atmospheric" [CO2*], air-sea CO2 flux, sfc DIC rate of change, & Delta pCO2
        co2starair = K0 * DBLE(fco2atm(i)) * 1.0e-6_r8 * DBLE(rhoSW(i)) !Equil. [CO2*] for atm CO2 at Patm & sfc-water T,S [mol/m3]
        co2star = DBLE(co2(i))                                          !Oceanic [CO2*] in [mol/m3] from vars.f90
        co2flux(i) = REAL(kwco2 * (co2starair - co2star))               !Air-sea CO2 flux [mol/(m2 * s)]
        co2ex(i) = co2flux(i) / dz1                                     !Change in sfc DIC due to gas exchange [mol/[m3 * s)]
        dpco2(i) = pco2(i) - pco2atm(i)                                 !Delta pCO2 (oceanic - atmospheric pCO2) [uatm]
     ENDIF

  END DO

  RETURN
END SUBROUTINE flxco2

!>    Compute xCO2 from arrays of pCO2atm, in situ T, S, & atm pressure
SUBROUTINE pCO2atm2xCO2(pCO2atm, temp, salt, Patm, N, xCO2)
  !    Purpose:
  !    Compute xCO2 from arrays of pCO2atm, in situ T, S, & atm pressure

  USE msingledouble

  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> atmospheric partial pressure of CO2 [uatm] 
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: pCO2atm
  !> in situ temperature [C]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: salt
  !> atmospheric pressure [atm]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: Patm
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> mole fraction of CO2 [ppm]
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: xCO2

! LOCAL variables:
  REAL(kind=r8) :: dpCO2atm, dPatm
  REAL(kind=r8), DIMENSION(N) :: pH20
  REAL(kind=r8) :: dxCO2

  INTEGER :: i

  call vapress(temp, salt, N, pH20)

  DO i = 1,N
     dpCO2atm  = DBLE(pCO2atm(i))
     dPatm     = DBLE(Patm(i))
     dxCO2     = dpCO2atm / (dPatm - pH20(i))
     xCO2(i) = REAL(dxCO2)
  END DO

  RETURN
END SUBROUTINE pCO2atm2xCO2

!>    Compute Schmidt number for CO2 in seawater from temperature
FUNCTION scco2(temp)

!  Computes the Schmidt number of CO2 in seawater using the
!  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
!  7373-7382).  Input is temperature in deg C.

   USE msingledouble
   IMPLICIT NONE

!  Input & output variables:
   REAL(kind=r8), INTENT(in) :: temp
!  REAL(kind=r8), INTENT(out) :: scco2
   REAL(kind=r8) :: scco2

   scco2 = 2073.1 - 125.62*temp + 3.6276*temp**2 - 0.043219*temp**3

  RETURN
END FUNCTION scco2

!>    Compute pCO2atm from arrays of xCO2, in situ T, S, & atm pressure
SUBROUTINE x2pCO2atm(xCO2, temp, salt, Patm, N, pCO2atm)
  !    Purpose:
  !    Compute pCO2atm from arrays of xCO2, in situ T, S, & atm pressure

  USE msingledouble

  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> mole fraction of CO2 [ppm]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: xCO2
  !> in situ temperature [C]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: salt
  !> atmospheric pressure [atm]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: Patm
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> oceanic partial pressure of CO2 [uatm] 
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: pCO2atm

! LOCAL variables:
  REAL(kind=r8) :: dxCO2, dPatm
  REAL(kind=r8), DIMENSION(N) :: pH20
  REAL(kind=r8) :: dpCO2atm

  INTEGER :: i

! Compute vapor pressure of seawater [in atm]
  call vapress(temp, salt, N, pH20)

  DO i = 1,N
     dxCO2     = DBLE(xCO2(i))
     dPatm     = DBLE(Patm(i))
     dpCO2atm = (dPatm - pH20(i)) * dxCO2
     pCO2atm(i) = REAL(dpCO2atm)
  END DO

  RETURN
END SUBROUTINE x2pCO2atm

!>    Compute vapor pressure of seawater (atm) following preocedure from Weiss & Price (1980)
SUBROUTINE vapress(temp, salt, N, vpsw)
  !    Purpose:
  !    Compute vapor pressure of seawater (atm) following preocedure from Weiss & Price (1980)

  USE msingledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> in situ temperature [C]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: salt
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> vapor pressure of seawater [atm] 
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: vpsw

! LOCAL variables:
  REAL(kind=r8) :: tk, dsalt

  INTEGER :: i

  DO i = 1,N
     dsalt = DBLE(salt(i))
     tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
     vpsw(i) = exp(24.4543d0 - 67.4509d0*(100.0d0/tk) - 4.8489d0*log(tk/100) - 0.000544d0*dsalt)
  END DO

  RETURN
END SUBROUTINE vapress

END MODULE gasx
