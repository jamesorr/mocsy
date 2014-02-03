!> Module with vars subroutine - compute carbonate system vars from DIC,Alk,T,S,P,nuts
MODULE mvars
CONTAINS
!>    Computes standard carbonate system variables (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
!!    as 1D arrays FROM
!!    temperature, salinity, pressure,
!!    total alkalinity (ALK), dissolved inorganic carbon (DIC),
!!    silica and phosphate concentrations (all 1-D arrays)
SUBROUTINE vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
                temp, sal, alk, dic, sil, phos, depth, lat, N,                          &
                optCON, optT, optP, optB, optK1K2, optKf                                  )

  !   Purpose:
  !     Computes other standard carbonate system variables (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
  !     as 1D arrays
  !     FROM:
  !     temperature, salinity, pressure,
  !     total alkalinity (ALK), dissolved inorganic carbon (DIC),
  !     silica and phosphate concentrations (all 1-D arrays)

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
  !     ----------
  !     optKf:
  !     ----------
  !       -> 'pf' means use Perez & Fraga (1987) formulation for Kf (recommended by Dickson et al., 2007)
  !               **** BUT Valid for  9 < T < 33°C and 10 < S < 40.
  !       -> 'dg' means use Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)

  !     OUTPUT variables:
  !     =================
  !     ph   = pH on total scale
  !     pco2 = CO2 partial pressure (uatm)
  !     fco2 = CO2 fugacity (uatm)
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
  USE mconstants
  USE mp80
  USE mrho
  USE msw_temp

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
  CHARACTER(3), INTENT(in) :: optB
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
  CHARACTER(2), INTENT(in) :: optKf
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) 
  CHARACTER(3), INTENT(in) :: optK1K2

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
  REAL(kind=r8) :: dfCO2, dpCO2
  REAL(kind=r8) :: drho
  REAL(kind=r8) :: invtk,dlogtk

  REAL(kind=r8) :: K0, K1, K2, Kb, Kw, Ks, Kf, Kspc
  REAL(kind=r8) :: Kspa, K1p, K2p, K3p, Ksi
  REAL(kind=r8) :: St, Ft, Bt

  REAL(kind=r8), DIMENSION(1) :: aK0, aK1, aK2, aKb, aKw, aKs, aKf, aKspc
  REAL(kind=r8), DIMENSION(1) :: aKspa, aK1p, aK2p, aK3p, aKsi
  REAL(kind=r8), DIMENSION(1) :: aSt, aFt, aBt

  REAL(kind=r8) :: DD,A,B,C,PhiD

  INTEGER :: i, icount

  REAL(kind=r8) :: t, tk, prb
  REAL(kind=r8) :: s
  REAL(kind=r8) :: tc, ta
  REAL(kind=r8) :: sit, pt
  REAL(kind=r8) :: ah1

  REAL(kind=r8) :: HSO4, HF, HSI, HPO4
  REAL(kind=r8) :: ab, aw, ac, ah2, erel

  REAL(kind=r8) :: cu, cb, cc
  INTEGER :: kcomp

! REAL(kind=r4) :: p80, sw_temp, rho
! EXTERNAL  p80, sw_temp, rho

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
     IF (temp(i) < -3. .OR. temp(i) > 100.) THEN
        print *,'OUT of RANGE TEMPERATURE: i, temp(i) =', i, temp(i)
     ENDIF
     IF (trim(optT) == 'Tpot') THEN
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
     ELSEIF (trim(optT) == 'Tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis(i) = temp(i)
        tempis68 = (temp(i) - 0.0002) / 0.99975
!       tempis(i) = 0.99975*tempis68 + 0.0002
     ELSE
        PRINT *,"optT must be either 'Tpot' or 'Tinsitu'"
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
        t = DBLE(tempis(i))
        tk = 273.15d0 + t
        invtk=1.0d0/tk
        dlogtk=LOG(tk)

!       Pressure effect (prb is in bars)
        prb = DBLE(p(i)) / 10.0d0

!       Salinity (equivalent array in double precision)
        s = DBLE(ssal)

!       Get all equilibrium constants and total concentrations of SO4, F, B
        CALL constants(aK0, aK1, aK2, aKb, aKw, aKs, aKf, aKspc, aKspa,  &
                       aK1p, aK2p, aK3p, aKsi,                           &
                       aSt, aFt, aBt,                                    &
                       temp(i), sal(i),                                  &
                       depth(i), lat(i), 1,                              &
                       optT, optP, optB, optK1K2, optKf)

!       Unlike f77, in F90 we can't assign an array (dimen=1) to a scalar in a routine argument
!       Thus, set scalar constants equal to array (dimension=1) values required as arguments
        K0 = aK0(1) ; K1 = aK1(1) ; K2 = aK2(1) ; Kb = aKb(1) ; Kw = aKw(1) 
        Ks = aKs(1) ; Kf = aKs(1) ; Kspc = aKspc(1) ; Kspa = aKspa(1) 
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

        tc = DBLE(sdic)/drho
        ta = DBLE(salk)/drho
        sit = DBLE(ssil)/drho
        pt  = DBLE(sphos)/drho

!       Computation loop for [H+]
!       Based on HAMOCC3 code (Maier-Reimer, 1993, Global Biogeoch. Cycles)
!       except (i) it is not made in a do loop,  
!             (ii) alk eqn includes contributions from P & Si (as in OCMIP2), &
!            (iii) the convergence criteria is slightly different
        ah1 = 1.e-8_r8
        kcomp = 0
2       kcomp = kcomp+1
!       print *, "kcomp  =  ", kcomp
        HSO4 = St/(1.0d0 + Ks/(ah1/(1.0d0 + St/Ks)))
        HF = 1.0d0/(1.0d0 + Kf/ah1)
        HSI = 1.0d0/(1.0d0 + ah1/Ksi)
        HPO4 = (K1p*K2p*(ah1 + 2.*K3p) - ah1**3)/   &
               (ah1**3 + K1p*ah1**2 + K1p*K2p*ah1 +   &
               K1p*K2p*K3p)

        ab = Bt/(1.0d0 + ah1/Kb)
        aw = Kw/ah1 - ah1/(1.0d0 + St/Ks)
        ac = ta + hso4 - sit*hsi - ab - aw + Ft*hf - pt*hpo4
        ah2 = SQRT((tc - ac)**2 + 4.0d0*(ac*K2/K1)*(2.0d0*tc - ac))
        ah2 = 0.5d0*K1/ac*((tc - ac) + ah2)
        erel = (ah2 - ah1)/ah2

        IF (ABS(erel) >= 1.e-7_r8) THEN
           ah1 = ah2
           IF (kcomp < 25) GOTO 2
        END IF

        IF (ah1 > 0.d0) THEN
           ph(i) = REAL(-LOG10(ah1))
        ELSE
           ph(i) = 1.e20_r4
        ENDIF

!       Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
        cu = (2.0d0 * tc - ac) / (2.0d0 + K1 / ah1)
        cb = K1 * cu / ah1
        cc = K2 * cb / ah1

!       If optCON = 'mol/m3', then:
!       convert output var concentrations from mol/kg to mol/m^3
!       e.g., for case when drho = 1028, multiply by [1.028 kg/L  x  1000 L/m^3])
        co2(i)  = REAL(cu * drho)
        hco3(i) = REAL(cb * drho)
        co3(i)  = REAL(cc * drho)

!       Determine CO2 pressure and fugacity (in microatm)
!       NOTE: equation below for pCO2 requires CO2* in mol/kg
        dfCO2 = cu   * 1.e6_r8/K0
        B=(-1636.75d0+12.0408d0*tk-0.0327957d0*(tk*tk)+0.0000316528d0* &
             (tk*tk*tk))*1e-6

        dpCO2= dfCO2 / EXP(100000.0d0 &
             *(B+ 2.0d0*(57.7d0 - 0.118d0*tk)*1e-6_r8) &
             / (8.314d0*tk) &
             )

        fCO2(i) = REAL(dfCO2)
        pCO2(i) = REAL(dpCO2)

!       Determine Omega Calcite et Aragonite
        OmegaA(i) = REAL(((0.01028d0*s/35.0d0)*cc)/Kspa)
        OmegaC(i) = REAL(((0.01028d0*s/35.0d0)*cc)/Kspc)

!       Determine Revelle Factor (BetaD w/ code from Seacarb software)
!       - analytical solution from Frankignoulle (1974)
        DD = -((-Kb*Bt)/((ah1 + Kb)*(ah1+Kb))) - (-Kw/(ah1*ah1)) + 1.0d0
        A = (2.0d0*K2*(2.0d0*cc + cb) + ah1*(ah1 + 2.0d0*K2)*DD)             &
            / ((ah1 + 2.0d0*K2)*(ah1 + 2.0d0*K2))
        B = ( ( (2.0d0*cc + cb) * ah1)/((ah1 + 2.0d0*K2)*K1) + (ah1/K1)* A )
        C = (-K2*(2.0d0*cc + cb) + K2*(2.0d0*K2 + ah1)*DD)                   &
            / ((ah1+2.0d0*K2)*(ah1+2.0d0*K2))
        PhiD = -1.0d0/(ah1*LOG(10.0d0) * ( B+A+C ) )
        BetaD(i) = REAL(-ah1*LOG(10.0d0)*tc/cu*B*PhiD)

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
