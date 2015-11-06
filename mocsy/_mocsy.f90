!> \file singledouble.f90
!! \BRIEF
!> Module that defines single and double precision - used by all other modules
MODULE msingledouble
! INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)
! INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12)
  INTEGER, PARAMETER :: r4 = KIND(1.0)
  INTEGER, PARAMETER :: r8 = KIND(1.0d0)
  INTEGER, PARAMETER :: wp = KIND(1.0d0)
END MODULE msingledouble


!> \file p80.f90
!! \BRIEF
!> Module with p80 function - compute pressure from depth
MODULE mp80
CONTAINS
!>     Function to compute pressure from depth using Saunder's (1981) formula with eos80.
FUNCTION p80(dpth,xlat)

  !     Compute Pressure from depth using Saunder's (1981) formula with eos80.

  !     Reference:
  !     Saunders, Peter M. (1981) Practical conversion of pressure
  !     to depth, J. Phys. Ooceanogr., 11, 573-574, (1981)

  !     Coded by:
  !     R. Millard
  !     March 9, 1983
  !     check value: p80=7500.004 dbars at lat=30 deg., depth=7321.45 meters

  !     Modified (slight format changes + added ref. details):
  !     J. Orr, 16 April 2009

  USE msingledouble
  IMPLICIT NONE

! Input variables:
  !> depth [m]
  REAL(kind=r4), INTENT(in) :: dpth
  !> latitude [degrees]
  REAL(kind=r4), INTENT(in) :: xlat

! Output variable:
  !> pressure [db]
  REAL(kind=r4) :: p80

! Local variables:
  REAL(kind=r4) :: pi
  REAL(kind=r4) :: plat, d, c1

  pi=3.141592654

  plat = ABS(xlat*pi/180.)
  d  = SIN(plat)
  c1 = 5.92e-3+d**2 * 5.25e-3

  p80 = ((1-c1)-SQRT(((1-c1)**2)-(8.84e-6*dpth))) / 4.42e-6

  RETURN
END FUNCTION p80
END MODULE mp80


!> \file sw_adtg.f90
!! \BRIEF
!> Module with sw_adtg function - compute adiabatic temp. gradient from S,T,P
MODULE msw_adtg
CONTAINS
!>  Function to calculate adiabatic temperature gradient as per UNESCO 1983 routines.
FUNCTION sw_adtg  (s,t,p)

  !     ==================================================================
  !     Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================
  USE msingledouble
  IMPLICIT NONE
  !> salinity [psu (PSU-78)]
  REAL(kind=r8) :: s
  !> temperature [degree C (IPTS-68)]
  REAL(kind=r8) :: t
  !> pressure [db]
  REAL(kind=r8) :: p

  REAL(kind=r8) :: a0,a1,a2,a3,b0,b1,c0,c1,c2,c3,d0,d1,e0,e1,e2
  REAL(kind=r8) :: sref

  REAL(kind=r8) :: sw_adtg

  sref = 35.d0
  a0 =  3.5803d-5
  a1 = +8.5258d-6
  a2 = -6.836d-8
  a3 =  6.6228d-10

  b0 = +1.8932d-6
  b1 = -4.2393d-8

  c0 = +1.8741d-8
  c1 = -6.7795d-10
  c2 = +8.733d-12
  c3 = -5.4481d-14

  d0 = -1.1351d-10
  d1 =  2.7759d-12

  e0 = -4.6206d-13
  e1 = +1.8676d-14
  e2 = -2.1687d-16

  sw_adtg =  a0 + (a1 + (a2 + a3*T)*T)*T &
       + (b0 + b1*T)*(S-sref) &
       + ( (c0 + (c1 + (c2 + c3*T)*T)*T) + (d0 + d1*T)*(S-sref) )*P &
       + (  e0 + (e1 + e2*T)*T )*P*P

END FUNCTION sw_adtg
END MODULE msw_adtg

!> \file sw_ptmp.f90
!! \BRIEF
!> Module with sw_ptmp function - compute potential T from in-situ T
MODULE msw_ptmp
CONTAINS
!> Function to calculate potential temperature [C] from in-situ temperature
FUNCTION sw_ptmp  (s,t,p,pr)

  !     ==================================================================
  !     Calculates potential temperature [C] from in-situ Temperature [C]
  !     From UNESCO 1983 report.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================

  !     Input arguments:
  !     -------------------------------------
  !     s  = salinity            [psu      (PSS-78) ]
  !     t  = temperature         [degree C (IPTS-68)]
  !     p  = pressure            [db]
  !     pr = reference pressure  [db]

  USE msingledouble
  USE msw_adtg
  IMPLICIT NONE

! Input arguments
  !> salinity [psu (PSS-78)]
  REAL(kind=r8) :: s
  !> temperature [degree C (IPTS-68)]
  REAL(kind=r8) :: t
  !> pressure [db]
  REAL(kind=r8) :: p
  !> reference pressure  [db]
  REAL(kind=r8) :: pr

! local arguments
  REAL(kind=r8) :: del_P ,del_th, th, q
  REAL(kind=r8) :: onehalf, two, three
  PARAMETER (onehalf = 0.5d0, two = 2.d0, three = 3.d0 )

! REAL(kind=r8) :: sw_adtg
! EXTERNAL sw_adtg

! Output
  REAL(kind=r8) :: sw_ptmp

  ! theta1
  del_P  = PR - P
  del_th = del_P*sw_adtg(S,T,P)
  th     = T + onehalf*del_th
  q      = del_th

  ! theta2
  del_th = del_P*sw_adtg(S,th,P+onehalf*del_P)
  th     = th + (1.d0 - 1.d0/SQRT(two))*(del_th - q)
  q      = (two-SQRT(two))*del_th + (-two+three/SQRT(two))*q

  ! theta3
  del_th = del_P*sw_adtg(S,th,P+onehalf*del_P)
  th     = th + (1.d0 + 1.d0/SQRT(two))*(del_th - q)
  q      = (two + SQRT(two))*del_th + (-two-three/SQRT(two))*q

  ! theta4
  del_th = del_P*sw_adtg(S,th,P+del_P)
  sw_ptmp     = th + (del_th - two*q)/(two*three)

  RETURN
END FUNCTION sw_ptmp
END MODULE msw_ptmp


!> \file sw_temp.f90
!! \BRIEF
!> Module with sw_temp function - compute in-situ T from potential T
MODULE msw_temp
CONTAINS
!> Function to compute in-situ temperature [C] from potential temperature [C]
FUNCTION sw_temp( s, t, p, pr )
  !     =============================================================
  !     SW_TEMP
  !     Computes in-situ temperature [C] from potential temperature [C]
  !     Routine available in seawater.f (used for MIT GCM)
  !     Downloaded seawater.f (on 17 April 2009) from
  !     http://ecco2.jpl.nasa.gov/data1/beaufort/MITgcm/bin/
  !     =============================================================

  !     REFERENCES:
  !     Fofonoff, P. and Millard, R.C. Jr
  !     Unesco 1983. Algorithms for computation of fundamental properties of
  !     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
  !     Eqn.(31) p.39

  !     Bryden, H. 1973.
  !     "New Polynomials for thermal expansion, adiabatic temperature gradient
  !     and potential temperature of sea water."
  !     DEEP-SEA RES., 1973, Vol20,401-408.
  !     =============================================================

  !     Simple modifications: J. C. Orr, 16 April 2009
  !     - combined fortran code from MITgcm site & simplification in
  !       CSIRO code (matlab equivalent) from Phil Morgan

  USE msingledouble
  USE msw_ptmp
  IMPLICIT NONE

  !     Input arguments:
  !     -----------------------------------------------
  !     s  = salinity              [psu      (PSS-78) ]
  !     t  = potential temperature [degree C (IPTS-68)]
  !     p  = pressure              [db]
  !     pr = reference pressure    [db]

  !> salinity [psu (PSS-78)]
  REAL(kind=r4) ::   s
  !> potential temperature [degree C (IPTS-68)]
  REAL(kind=r4) ::   t
  !> pressure [db]
  REAL(kind=r4) ::   p
  !> reference pressure [db]
  REAL(kind=r4) ::   pr

  REAL(kind=r8) ::  ds, dt, dp, dpr
  REAL(kind=r8) :: dsw_temp

  REAL(kind=r4) ::   sw_temp
! EXTERNAL sw_ptmp
! REAL(kind=r8) ::   sw_ptmp

  ds = DBLE(s)
  dt = DBLE(t)
  dp = DBLE(p)
  dpr = DBLE(pr)

  !    Simple solution
  !    (see https://svn.mpl.ird.fr/us191/oceano/tags/V0/lib/matlab/seawater/sw_temp.m)
  !    Carry out inverse calculation by swapping P_ref (pr) and Pressure (p)
  !    in routine that is normally used to compute potential temp from temp
  dsw_temp = sw_ptmp(ds, dt, dpr, dp)
  sw_temp = REAL(dsw_temp)

  !    The above simplification works extremely well (compared to Table in 1983 report)
  !    whereas the sw_temp routine from MIT GCM site does not seem to work right

  RETURN
END FUNCTION sw_temp
END MODULE msw_temp


!> \file constants.f90
!! \BRIEF
!> Module with contants subroutine - computes carbonate system constants
!! from S,T,P
MODULE mconstants
CONTAINS
!> Compute thermodynamic constants
!! FROM temperature, salinity, and pressure (1D arrays)
SUBROUTINE constants(K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa,  &
                     K1p, K2p, K3p, Ksi,                      &
                     St, Ft, Bt,                              &
                     temp, sal, Patm,                         &
                     depth, lat, N,                           &
                     optT, optP, optB, optK1K2, optKf, optGAS)

  !   Purpose:
  !     Compute thermodynamic constants
  !     FROM: temperature, salinity, and pressure (1D arrays)

  !     INPUT variables:
  !     ================
  !     Patm    = atmospheric pressure [atm]
  !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP='db')
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
  !             = dummy array (unused when optP='db')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !     sal     = salinity in [psu]
  !     ---------
  !     optT: choose in situ vs. potential temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
  !       -> 'Tinsitu' means input is already in-situ Temperature, not pot. Temp ("tempis" not computed)
  !     ---------
  !     optP: choose depth (m) vs pressure (db) as input
  !     ---------
  !       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
  !       -> 'db' means "depth" input is already in situ pressure [db], not m (p = depth)
  !     ---------
  !     optB:
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
  !     -----------
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
  !     K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi
  !     St, Ft, Bt

  USE msingledouble
  USE mp80
  USE msw_temp
  USE msw_ptmp
  IMPLICIT NONE

! Input variables
  !>     number of records
  INTEGER, INTENT(in) :: N
  !> in <b>situ temperature</b> (when optT='Tinsitu', typical data)
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) [degree C]
  REAL(kind=r4), INTENT(in),    DIMENSION(N) :: temp
  !> depth in <b>meters</b> (when optP='m') or <b>decibars</b> (when optP='db')
  REAL(kind=r4), INTENT(in),    DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=r4), INTENT(in),    DIMENSION(N) :: lat
  !> salinity <b>[psu]</b>
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: sal
!f2py optional , depend(sal) :: n=len(sal)

  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: Patm

  !> for temp input, choose \b 'Tinsitu' for in situ Temp or
  !! \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
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

! Ouput variables
  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K0
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K2
  !> equilibrium constant for dissociation of boric acid
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride
  !! either from Dickson and Riley (1979) or from Perez and Fraga (1987), depending on optKf
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Ksi
  !> total sulfate (Morris & Riley, 1966)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: St
  !> total fluoride  (Riley, 1965)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Ft
  !> total boron
  !! from either Uppstrom (1974) or Lee et al. (2010), depending on optB
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Bt

! Local variables
  REAL(kind=r4) :: ssal
  REAL(kind=r4) :: p
  REAL(kind=r4) :: tempot, tempis68, tempot68
  REAL(kind=r4) :: tempis
  REAL(kind=r8) :: is, invtk, dlogtk, is2, s2, sqrtis
  REAL(kind=r8) :: Ks_0p, Kf_0p
  REAL(kind=r8) :: total2free, free2SWS, total2SWS, SWS2total
  REAL(kind=r8) :: total2free_0p, free2SWS_0p, total2SWS_0p
! REAL(kind=r8) :: free2SWS, free2SWS_0p

  REAL(kind=r8) :: dtempot, dtempot68
  REAL(kind=r8) :: R

  REAL(kind=r8) :: pK1o, ma1, mb1, mc1, pK1
  REAL(kind=r8) :: pK2o, ma2, mb2, mc2, pK2

  REAL(kind=r8), DIMENSION(12) :: a0, a1, a2, b0, b1, b2
  REAL(kind=r8), DIMENSION(12) :: deltav, deltak, lnkpok0
  REAL(kind=r8) :: tmp, nK0we74

  INTEGER :: i, icount, ipc

  REAL(kind=r8) :: t, tk, tk0, prb
  REAL(kind=r8) :: s, sqrts, s15, scl

  REAL(kind=r8) :: Phydro_atm, Patmd, Ptot, Rgas_atm, vbarCO2

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS

  ! CONSTANTS
  ! =========
  ! Constants in formulation for Pressure effect on K's (Millero, 95)
  ! with corrected coefficients for Kb, Kw, Ksi, etc.

  ! index: 1) K1 , 2) K2, 3) Kb, 4) Kw, 5) Ks, 6) Kf, 7) Kspc, 8) Kspa,
  !            9) K1P, 10) K2P, 11) K3P, 12) Ksi

  DATA a0 /-25.5_r8, -15.82_r8, -29.48_r8, -20.02_r8, &
          -18.03_r8,  -9.78_r8, -48.76_r8, -45.96_r8, &
          -14.51_r8, -23.12_r8, -26.57_r8, -29.48_r8/
  DATA a1 /0.1271_r8, -0.0219_r8, 0.1622_r8, 0.1119_r8, &
           0.0466_r8, -0.0090_r8, 0.5304_r8, 0.5304_r8, &
           0.1211_r8, 0.1758_r8, 0.2020_r8, 0.1622_r8/
  DATA a2 /     0.0_r8,       0.0_r8, -2.608e-3_r8, -1.409e-3_r8, &
           0.316e-3_r8, -0.942e-3_r8,  0.0_r8,       0.0_r8, &
          -0.321e-3_r8, -2.647e-3_r8, -3.042e-3_r8, -2.6080e-3_r8/
  DATA b0 /-3.08e-3_r8, 1.13e-3_r8,  -2.84e-3_r8,   -5.13e-3_r8, &
           -4.53e-3_r8, -3.91e-3_r8, -11.76e-3_r8, -11.76e-3_r8, &
           -2.67e-3_r8, -5.15e-3_r8,  -4.08e-3_r8,  -2.84e-3_r8/
  DATA b1 /0.0877e-3_r8, -0.1475e-3_r8, 0.0_r8,       0.0794e-3_r8, &
           0.09e-3_r8,    0.054e-3_r8,  0.3692e-3_r8, 0.3692e-3_r8, &
           0.0427e-3_r8,  0.09e-3_r8,   0.0714e-3_r8, 0.0_r8/
  DATA b2 /12*0.0_r8/

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with
!        the !f2py statements that precede each type declaraion
  IF (PRESENT(optB)) THEN
    opB = optB
  ELSE
    opB = 'l10'
  ENDIF
  IF (PRESENT(optKf)) THEN
    opKf = optKf
  ELSE
    opKf = 'pf'
  ENDIF
  IF (PRESENT(optK1K2)) THEN
    opK1K2 = optK1K2
  ELSE
    opK1K2 = 'l'
  ENDIF
  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF

  R = 83.14472_r8

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
!     - depth [m] => convert to pressure [db]
!     - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
     IF (trim(optP) == 'm' ) THEN
!       Compute pressure [db] from depth [m] and latitude [degrees]
        p = p80(depth(i), lat(i))
     ELSEIF (trim(optP) == 'db' ) THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        p = depth(i)
     ELSE
        PRINT *,"optP must be 'm' or 'db'"
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
        tempis68 = sw_temp(sal(i), tempot68, p, 0. )
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis = 0.99975*tempis68 + 0.0002
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (trim(optT) == 'Tinsitu' .OR. trim(optT) == 'tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis    = temp(i)
        tempis68  = (temp(i) - 0.0002) / 0.99975
        dtempot68 = sw_ptmp(DBLE(sal(i)), DBLE(tempis68), DBLE(p), 0.0d0)
        dtempot   = 0.99975*dtempot68 + 0.0002
     ELSE
        PRINT *,"optT must be either 'Tpot' or 'Tinsitu'"
        PRINT *,"you specified optT =", trim(optT)
        STOP
     ENDIF

!    Compute constants:
     IF (temp(i) >= -5. .AND. temp(i) < 1.0e+2) THEN
!       Test to indicate if any of input variables are unreasonable
        IF (      sal(i) < 0.  .OR.  sal(i) > 1e+3) THEN
           PRINT *, 'i, icount, temp, sal =', i, icount, temp(i), sal(i)
        ENDIF
!       Zero out negative salinity (prev case for OCMIP2 model w/ slightly negative S in some coastal cells)
        IF (sal(i) < 0.0) THEN
           ssal = 0.0
        ELSE
           ssal = sal(i)
        ENDIF

!       Absolute temperature (Kelvin) and related values
        t = DBLE(tempis)
        tk = 273.15d0 + t
        invtk=1.0d0/tk
        dlogtk=LOG(tk)

!       Atmospheric pressure
        Patmd = DBLE(Patm(i))

!       Hydrostatic pressure (prb is in bars)
        prb = DBLE(p) / 10.0d0

!       Salinity and simply related values
        s = DBLE(ssal)
        s2=s*s
        sqrts=SQRT(s)
        s15=s**1.5d0
        scl=s/1.80655d0

!       Ionic strength:
        is = 19.924d0*s/(1000.0d0 - 1.005d0*s)
        is2 = is*is
        sqrtis = SQRT(is)

!       Total concentrations for sulfate, fluoride, and boron

!       Sulfate: Morris & Riley (1966)
        St(i) = 0.14d0 * scl/96.062d0

!       Fluoride:  Riley (1965)
        Ft(i) = 0.000067d0 * scl/18.9984d0

!       Boron:
        IF (trim(opB) == 'l10') THEN
!          New formulation from Lee et al (2010)
           Bt(i) = 0.0002414d0 * scl/10.811d0
        ELSEIF (trim(opB) == 'u74') THEN
!          Classic formulation from Uppström (1974)
           Bt(i) = 0.000232d0  * scl/10.811d0
        ELSE
           PRINT *,"optB must be 'l10' or 'u74'"
           STOP
        ENDIF

!       K0 (K Henry)
!       CO2(g) <-> CO2(aq.)
!       K0  = [CO2]/ fCO2
!       Weiss (1974)   [mol/kg/atm]
        IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
           tk0 = tk                   !in situ temperature (K) for K0 calculation
           Ptot = Patmd               !total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Ppot'    .OR. trim(opGAS) == 'ppot') THEN
           tk0 = dtempot + 273.15d0   !potential temperature (K) for K0 calculation as needed for potential fCO2 & pCO2
           Ptot = Patmd               !total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
           tk0 = tk                     !in situ temperature (K) for K0 calculation
           Phydro_atm = prb / 1.01325d0 !convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
           Ptot = Patmd + Phydro_atm    !total pressure (in atm) = atmospheric pressure + hydrostatic pressure
        ELSE
           PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
           STOP
        ENDIF
        tmp = 9345.17d0/tk0 - 60.2409d0 + 23.3585d0 * LOG(tk0/100.0d0)
        nK0we74 = tmp + s*(0.023517d0 - 0.00023656d0*tk0 + 0.0047036e-4_r8*tk0*tk0)
        K0(i) = EXP(nK0we74)

!       K1 = [H][HCO3]/[H2CO3]
!       K2 = [H][CO3]/[HCO3]
        IF (trim(opK1K2) == 'l') THEN
!         Mehrbach et al. (1973) refit, by Lueker et al. (2000) (total pH scale)
          K1(i) = 10.0d0**(-1.0d0*(3633.86d0*invtk - 61.2172d0 + 9.6777d0*dlogtk  &
                  - 0.011555d0*s + 0.0001152d0*s2))
          K2(i) = 10.0d0**(-1*(471.78d0*invtk + 25.9290d0 - 3.16967d0*dlogtk      &
                  - 0.01781d0*s + 0.0001122d0*s2))
        ELSEIF (trim(opK1K2) == 'm10') THEN
!         Millero (2010, Mar. Fresh Wat. Res.) (total pH scale)
!         pK1o = 6320.813d0*invtk + 19.568224d0*dlogtk -126.34048d0
!         ma1 = 13.4051d0*sqrts + 0.03185d0*s - (5.218e-5)*s2
!         mb1 = -531.095d0*sqrts - 5.7789d0*s
!         mc1 = -2.0663d0*sqrts
!         pK1 = pK1o + ma1 + mb1*invtk + mc1*dlogtk
!         K1(i) = 10.0d0**(-pK1)

!         pK2o = 5143.692d0*invtk + 14.613358d0*dlogtk -90.18333d0
!         ma2 = 21.5724d0*sqrts + 0.1212d0*s - (3.714e-4)*s2
!         mb2 = -798.292d0*sqrts - 18.951d0*s
!         mc2 = -3.403d0*sqrts
!         pK2 = pK2o + ma2 + mb2*invtk + mc2*dlogtk
!         K2(i) = 10.0d0**(-pK2)

!         Millero (2010, Mar. Fresh Wat. Res.) (seawater pH scale)
          pK1o = 6320.813d0*invtk + 19.568224d0*dlogtk -126.34048d0
          ma1 = 13.4038d0*sqrts + 0.03206d0*s - (5.242e-5)*s2
          mb1 = -530.659d0*sqrts - 5.8210d0*s
          mc1 = -2.0664d0*sqrts
          pK1 = pK1o + ma1 + mb1*invtk + mc1*dlogtk
          K1(i) = 10.0d0**(-pK1)

          pK2o = 5143.692d0*invtk + 14.613358d0*dlogtk -90.18333d0
          ma2 = 21.3728d0*sqrts + 0.1218d0*s - (3.688e-4)*s2
          mb2 = -788.289d0*sqrts - 19.189d0*s
          mc2 = -3.374d0*sqrts
          pK2 = pK2o + ma2 + mb2*invtk + mc2*dlogtk
          K2(i) = 10.0d0**(-pK2)
        ELSE
           PRINT *, "optK1K2 must be either 'l' or 'm10'"
           STOP
        ENDIF

!       Kb = [H][BO2]/[HBO2]
!       (total scale)
!       Millero p.669 (1995) using data from Dickson (1990)
        Kb(i) = EXP((-8966.90d0 - 2890.53d0*sqrts - 77.942d0*s +  &
                1.728d0*s15 - 0.0996d0*s2)*invtk +              &
                (148.0248d0 + 137.1942d0*sqrts + 1.62142d0*s) +   &
                (-24.4344d0 - 25.085d0*sqrts - 0.2474d0*s) *      &
                dlogtk + 0.053105d0*sqrts*tk)

!       K1p = [H][H2PO4]/[H3PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 65
!       Use Millero equation's 115.540 constant instead of 115.525 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K1p(i) = EXP(-4576.752d0*invtk + 115.540d0 - 18.453d0*dlogtk +  &
                 (-106.736d0*invtk + 0.69171d0) * sqrts +             &
                 (-0.65643d0*invtk - 0.01844d0) * s)

!       K2p = [H][HPO4]/[H2PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
!       Millero (1995), p.670, eq. 66
!       Use Millero equation's 172.1033 constant instead of 172.0833 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K2p(i) = EXP(-8814.715d0*invtk + 172.1033d0 - 27.927d0*dlogtk +  &
                 (-160.340d0*invtk + 1.3566d0)*sqrts +                 &
                 (0.37335d0*invtk - 0.05778d0)*s)

!       K3p = [H][PO4]/[HPO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 67
!       Use Millero equation's 18.126 constant instead of 18.141 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K3p(i) = EXP(-3070.75d0*invtk - 18.126d0 +            &
                 (17.27039d0*invtk + 2.81197d0) *             &
                 sqrts + (-44.99486d0*invtk - 0.09984d0) * s)

!       Ksi = [H][SiO(OH)3]/[Si(OH)4]
!       (seawater scale)
!       Millero (1995), p.671, eq. 72
!       Use Millero equation's 117.400 constant instead of 117.385 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        Ksi(i) = EXP(-8904.2d0*invtk  + 117.400d0 - 19.334d0*dlogtk +  &
                 (-458.79d0*invtk + 3.5913d0) * sqrtis +             &
                 (188.74d0*invtk - 1.5998d0) * is +                  &
                 (-12.1652d0*invtk + 0.07871d0) * is2 +              &
                 LOG(1.0 - 0.001005d0*s))

!       Kw = [H][OH]
!       (seawater scale)
!       Millero (1995) p.670, eq. 63 from composite data
!       Use Millero equation's 148.9802 constant instead of 148.9652 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        Kw(i) = EXP(-13847.26d0*invtk + 148.9802d0 - 23.6521d0*dlogtk +  &
               (118.67d0*invtk - 5.977d0 + 1.0495d0 * dlogtk) *          &
               sqrts - 0.01615d0 * s)

!       Ks = [H][SO4]/[HSO4]
!       (free scale)
!       Dickson (1990, J. chem. Thermodynamics 22, 113)
        Ks_0p = EXP(-4276.1d0*invtk + 141.328d0 - 23.093d0*dlogtk          &
                + (-13856.d0*invtk + 324.57d0 - 47.986d0*dlogtk) * sqrtis  &
                + (35474.d0*invtk - 771.54 + 114.723d0*dlogtk) * is      &
                - 2698.d0*invtk*is**1.5 + 1776.d0*invtk*is2              &
                + LOG(1.0d0 - 0.001005d0*s))

!       Kf = [H][F]/[HF]
!       (total scale)
        IF (trim(opKf) == 'dg') THEN
!          Dickson and Riley (1979) -- change pH scale to total (following Dickson & Goyet, 1994)
           Kf_0p = EXP(1590.2d0*invtk - 12.641d0 + 1.525d0*sqrtis +  &
                   LOG(1.0d0 - 0.001005d0*s) +                     &
                   LOG(1.0d0 + St(i)/Ks_0p))
        ELSEIF (trim(opKf) == 'pf') THEN
!          Perez and Fraga (1987) - Already on Total scale (no need for last line above)
!          Formulation as given in Dickson et al. (2007)
           Kf_0p = EXP(874.d0*invtk - 9.68d0 + 0.111d0*sqrts)
        ELSE
           PRINT *, "optKf must be either 'dg' or 'pf'"
           STOP
        ENDIF

!       Kspc (calcite) - apparent solubility product of calcite
!       (no scale)
!       Kspc = [Ca2+] [CO32-] when soln is in equilibrium w/ calcite
!       Mucci 1983 mol/kg-soln
        Kspc(i) = 10d0**(-171.9065d0 - 0.077993d0*tk + 2839.319d0/tk    &
                 + 71.595d0*LOG10(tk)                             &
                 + (-0.77712d0 + 0.0028426d0*tk + 178.34d0/tk)*sqrts  &
                 -0.07711d0*s + 0.0041249d0*s15 )


!       Kspa (aragonite) - apparent solubility product of aragonite
!       (no scale)
!       Kspa = [Ca2+] [CO32-] when soln is in equilibrium w/ aragonite
!       Mucci 1983 mol/kg-soln
        Kspa(i) = 10.d0**(-171.945d0 - 0.077993d0*tk + 2903.293d0/tk &
             +71.595d0*LOG10(tk) &
             +(-0.068393d0 + 0.0017276d0*tk + 88.135d0/tk)*sqrts &
             -0.10018d0*s + 0.0059415d0*s15 )

!       Pressure effect on K0 based on Weiss (1974, equation 5)
        Rgas_atm = 82.05736_r8      ! (cm3 * atm) / (mol * K)  CODATA (2006)
        vbarCO2 = 32.3_r8           ! partial molal volume (cm3 / mol) from Weiss (1974, Appendix, paragraph 3)
        K0(i) = K0(i) * exp( ((1-Ptot)*vbarCO2)/(Rgas_atm*tk0) )   ! Weiss (1974, equation 5)

!       Pressure effect on all other K's (based on Millero, (1995)
!           index: K1(1), K2(2), Kb(3), Kw(4), Ks(5), Kf(6), Kspc(7), Kspa(8),
!                  K1p(9), K2p(10), K3p(11), Ksi(12)
        DO ipc = 1, 12
           deltav(ipc)  =  a0(ipc) + a1(ipc) *t + a2(ipc) *t*t
           deltak(ipc)   = (b0(ipc)  + b1(ipc) *t + b2(ipc) *t*t)
           lnkpok0(ipc)  = (-(deltav(ipc)) &
                +(0.5d0*deltak(ipc) * prb) &
                )                         * prb/(R*tk)
        END DO

!       Pressure correction on Ks (Free scale)
        Ks(i) = Ks_0p*EXP(lnkpok0(5))
!       Conversion factor total -> free scale
        total2free     = 1.d0/(1.d0 + St(i)/Ks(i))   ! Kfree = Ktotal*total2free
!       Conversion factor total -> free scale at pressure zero
        total2free_0p  = 1.d0/(1.d0 + St(i)/Ks_0p)   ! Kfree = Ktotal*total2free

!       Pressure correction on Kf
!       Kf must be on FREE scale before correction
        Kf_0p = Kf_0p * total2free_0p   !Convert from Total to Free scale (pressure 0)
        Kf(i) = Kf_0p * EXP(lnkpok0(6)) !Pressure correction (on Free scale)
        Kf(i) = Kf(i)/total2free        !Convert back from Free to Total scale

!       Convert between seawater and total hydrogen (pH) scales
        free2SWS  = 1.d0 + St(i)/Ks(i) + Ft(i)/(Kf(i)*total2free)  ! using Kf on free scale
        total2SWS = total2free * free2SWS                          ! KSWS = Ktotal*total2SWS
        SWS2total = 1.d0 / total2SWS
!       Conversion at pressure zero
        free2SWS_0p  = 1.d0 + St(i)/Ks_0p + Ft(i)/(Kf_0p)  ! using Kf on free scale
        total2SWS_0p = total2free_0p * free2SWS_0p         ! KSWS = Ktotal*total2SWS

!       Convert from Total to Seawater scale before pressure correction
!       Must change to SEAWATER scale: K1, K2, Kb
        IF (trim(optK1K2) == 'l') THEN
          K1(i)  = K1(i)*total2SWS_0p
          K2(i)  = K2(i)*total2SWS_0p
          !This conversion is unnecessary for the K1,K2 from Millero (2010),
          !since we use here the formulation already on the seawater scale
        ENDIF
        Kb(i)  = Kb(i)*total2SWS_0p

!       Already on SEAWATER scale: K1p, K2p, K3p, Kb, Ksi, Kw

!       Other contants (keep on another scale):
!          - K0         (independent of pH scale, already pressure corrected)
!          - Ks         (already on Free scale;   already pressure corrected)
!          - Kf         (already on Total scale;  already pressure corrected)
!          - Kspc, Kspa (independent of pH scale; pressure-corrected below)

!       Perform actual pressure correction (on seawater scale)
        K1(i)   = K1(i)*EXP(lnkpok0(1))
        K2(i)   = K2(i)*EXP(lnkpok0(2))
        Kb(i)   = Kb(i)*EXP(lnkpok0(3))
        Kw(i)   = Kw(i)*EXP(lnkpok0(4))
        Kspc(i) = Kspc(i)*EXP(lnkpok0(7))
        Kspa(i) = Kspa(i)*EXP(lnkpok0(8))
        K1p(i)  = K1p(i)*EXP(lnkpok0(9))
        K2p(i)  = K2p(i)*EXP(lnkpok0(10))
        K3p(i)  = K3p(i)*EXP(lnkpok0(11))
        Ksi(i)  = Ksi(i)*EXP(lnkpok0(12))

!       Convert back to original total scale:
        K1(i)  = K1(i) *SWS2total
        K2(i)  = K2(i) *SWS2total
        K1p(i) = K1p(i)*SWS2total
        K2p(i) = K2p(i)*SWS2total
        K3p(i) = K3p(i)*SWS2total
        Kb(i)  = Kb(i) *SWS2total
        Ksi(i) = Ksi(i)*SWS2total
        Kw(i)  = Kw(i) *SWS2total

     ELSE

        K0(i)   = 1.e20_r8
        K1(i)   = 1.e20_r8
        K2(i)   = 1.e20_r8
        Kb(i)   = 1.e20_r8
        Kw(i)   = 1.e20_r8
        Ks(i)   = 1.e20_r8
        Kf(i)   = 1.e20_r8
        Kspc(i) = 1.e20_r8
        Kspa(i) = 1.e20_r8
        K1p(i)  = 1.e20_r8
        K2p(i)  = 1.e20_r8
        K3p(i)  = 1.e20_r8
        Ksi(i)  = 1.e20_r8
        Bt(i)   = 1.e20_r8
        Ft(i)   = 1.e20_r8
        St(i)   = 1.e20_r8

     ENDIF

  END DO

  RETURN
END SUBROUTINE constants
END MODULE mconstants


!> \file p2fCO2.f90
!! \BRIEF
!>    Module with p2fCO2 subroutine - compute fCO2 from pCO2, in situ T, atm pressure, hydrostatic pressure
MODULE mp2fCO2
CONTAINS
!>    Compute fCO2 from arrays of pCO2, in situ temp, atm pressure, & hydrostatic pressure
SUBROUTINE p2fCO2(pCO2, temp, Patm, p, N, fCO2)
  !    Purpose:
  !    Compute fCO2 from arrays of pCO2, in situ temp, atm pressure, & hydrostatic pressure

  USE msingledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, INTENT(in) :: N

! INPUT variables
  !> oceanic partial pressure of CO2 [uatm]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: pCO2
  !> in situ temperature [C]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: temp
  !> atmospheric pressure [atm]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: Patm
  !> hydrostatic pressure [db]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: p
!f2py optional , depend(pCO2) :: n=len(pCO2)

! OUTPUT variables:
  !> fugacity of CO2 [uatm]
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: fCO2

! LOCAL variables:
  REAL(kind=r8) :: dpCO2, dtemp, tk, dPatm, prb
  REAL(kind=r8) :: Ptot, Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  REAL(kind=r8) :: dfCO2

  INTEGER :: i

! REAL(kind=r8) :: sw_ptmp
! EXTERNAL sw_ptmp

  DO i = 1,N
     dpCO2     = DBLE(pCO2(i))
     dtemp     = DBLE(temp(i))
     dPatm     = DBLE(Patm(i))
     tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
     prb = DBLE(p(i)) / 10.0d0         !Pressure effect (prb is in bars)
     Ptot = dPatm + prb/1.01325d0      !Total pressure (atmospheric + hydrostatic) [atm]
     Rgas_atm = 82.05736_r8            !R in (cm3 * atm) / (mol * K)  from CODATA (2006)
!    To compute fugcoeff, we need 3 other terms (B, Del, xc2) as well as 3 others above (tk, Ptot, Rgas_atm)
     B = -1636.75d0 + 12.0408d0*tk - 0.0327957d0*(tk*tk) + 0.0000316528d0*(tk*tk*tk)
     Del = 57.7d0 - 0.118d0*tk
!    "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
!    x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
!    Let's assume that xCO2 = pCO2. Resulting fugcoeff is identical to 8th digit after the decimal.
     xCO2approx = dpCO2 * 1.e-6_r8
     xc2 = (1.0d0 - xCO2approx)**2
     fugcoeff = EXP( Ptot*(B + 2.0d0*xc2*Del)/(Rgas_atm*tk) )
     dfCO2 = dpCO2 * fugcoeff
     fCO2(i) = REAL(dfCO2)
  END DO

  RETURN
END SUBROUTINE p2fCO2
END MODULE mp2fCO2


!> \file phsolvers.f90
!! \BRIEF
!> Module with routines needed to solve pH-total alkalinity equation (Munhoven, 2013, GMD)
MODULE mphsolvers
!   Module of fastest solvers from Munhoven (2013, Geosci. Model Dev., 6, 1367-1388)
!   ! Taken from SolveSAPHE (mod_phsolvers.F90) & adapted very slightly for use with mocsy
!   ! SolveSaphe is distributed under the GNU Lesser General Public License
!   ! mocsy is distributed under the MIT License
!
! Modifications J. C. Orr, LSCE/IPSL, CEA-CNRS-UVSQ, France, 11 Sep 2014:
! 1) kept only the 3 fastest solvers (atgen, atsec, atfast) and routines which they call
! 2) reduced vertical white space: deleted many blank lines & comment lines that served as divisions
! 3) converted name from .F90 to .f90, deleting a few optional preprocesse if statements
! 4) read in mocsy computed equilibrium constants (as arguments) instead of USE MOD_CHEMCONST
! 5) converted routine names from upper case to lower case
! 6) commented out arguments and equations for NH4 and H2S acid systems

USE msingledouble
IMPLICIT NONE

! General parameters
REAL(KIND=wp), PARAMETER :: pp_rdel_ah_target = 1.E-8_wp
REAL(KIND=wp), PARAMETER :: pp_ln10 = 2.302585092994045684018_wp

! Maximum number of iterations for each method
INTEGER, PARAMETER :: jp_maxniter_atgen    = 50
INTEGER, PARAMETER :: jp_maxniter_atsec    = 50
INTEGER, PARAMETER :: jp_maxniter_atfast   = 50

! Bookkeeping variables for each method
! - SOLVE_AT_GENERAL
INTEGER :: niter_atgen    = jp_maxniter_atgen

! - SOLVE_AT_GENERAL_SEC
INTEGER :: niter_atsec    = jp_maxniter_atsec

! - SOLVE_AT_FAST (variant of SOLVE_AT_GENERAL w/o bracketing
INTEGER :: niter_atfast   = jp_maxniter_atfast

! Keep the following functions private to avoid conflicts with
! other modules that provide similar ones.
!PRIVATE AHINI_FOR_AT

CONTAINS
!===============================================================================
SUBROUTINE anw_infsup(p_dictot, p_bortot,                                     &
                      p_po4tot, p_siltot,                                     &
                      p_so4tot, p_flutot,                                     &
                      p_alknw_inf, p_alknw_sup)

! Subroutine returns the lower and upper bounds of "non-water-selfionization"
! contributions to total alkalinity (the infimum and the supremum), i.e
! inf(TA - [OH-] + [H+]) and sup(TA - [OH-] + [H+])

USE msingledouble
IMPLICIT NONE

! Argument variables
REAL(KIND=wp), INTENT(IN)  :: p_dictot
REAL(KIND=wp), INTENT(IN)  :: p_bortot
REAL(KIND=wp), INTENT(IN)  :: p_po4tot
REAL(KIND=wp), INTENT(IN)  :: p_siltot
!REAL(KIND=wp), INTENT(IN)  :: p_nh4tot
!REAL(KIND=wp), INTENT(IN)  :: p_h2stot
REAL(KIND=wp), INTENT(IN)  :: p_so4tot
REAL(KIND=wp), INTENT(IN)  :: p_flutot
REAL(KIND=wp), INTENT(OUT) :: p_alknw_inf
REAL(KIND=wp), INTENT(OUT) :: p_alknw_sup

p_alknw_inf =  -p_po4tot - p_so4tot - p_flutot
p_alknw_sup =   p_dictot + p_dictot + p_bortot &
              + p_po4tot + p_po4tot + p_siltot !&
!             + p_nh4tot + p_h2stot

RETURN
END SUBROUTINE anw_infsup

!===============================================================================

FUNCTION equation_at(p_alktot, p_h,       p_dictot, p_bortot,                 &
                     p_po4tot, p_siltot,                                      &
                     p_so4tot, p_flutot,                                      &
                     K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,          &
                     p_deriveqn)

USE msingledouble
IMPLICIT NONE
REAL(KIND=wp) :: equation_at

! Argument variables
REAL(KIND=wp), INTENT(IN)            :: p_alktot
REAL(KIND=wp), INTENT(IN)            :: p_h
REAL(KIND=wp), INTENT(IN)            :: p_dictot
REAL(KIND=wp), INTENT(IN)            :: p_bortot
REAL(KIND=wp), INTENT(IN)            :: p_po4tot
REAL(KIND=wp), INTENT(IN)            :: p_siltot
!REAL(KIND=wp), INTENT(IN)            :: p_nh4tot
!REAL(KIND=wp), INTENT(IN)            :: p_h2stot
REAL(KIND=wp), INTENT(IN)            :: p_so4tot
REAL(KIND=wp), INTENT(IN)            :: p_flutot
REAL(KIND=wp), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
REAL(KIND=wp), INTENT(IN)            :: K1p, K2p, K3p, Ksi
REAL(KIND=wp), INTENT(OUT), OPTIONAL :: p_deriveqn

! Local variables
!-----------------
REAL(KIND=wp) :: znumer_dic, zdnumer_dic, zdenom_dic, zalk_dic, zdalk_dic
REAL(KIND=wp) :: znumer_bor, zdnumer_bor, zdenom_bor, zalk_bor, zdalk_bor
REAL(KIND=wp) :: znumer_po4, zdnumer_po4, zdenom_po4, zalk_po4, zdalk_po4
REAL(KIND=wp) :: znumer_sil, zdnumer_sil, zdenom_sil, zalk_sil, zdalk_sil
REAL(KIND=wp) :: znumer_nh4, zdnumer_nh4, zdenom_nh4, zalk_nh4, zdalk_nh4
REAL(KIND=wp) :: znumer_h2s, zdnumer_h2s, zdenom_h2s, zalk_h2s, zdalk_h2s
REAL(KIND=wp) :: znumer_so4, zdnumer_so4, zdenom_so4, zalk_so4, zdalk_so4
REAL(KIND=wp) :: znumer_flu, zdnumer_flu, zdenom_flu, zalk_flu, zdalk_flu
REAL(KIND=wp) ::                                      zalk_wat, zdalk_wat
REAL(KIND=wp) :: aphscale

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
aphscale = 1._wp + p_so4tot/Ks

! H2CO3 - HCO3 - CO3 : n=2, m=0
znumer_dic = 2._wp*K1*K2 + p_h*       K1
zdenom_dic =       K1*K2 + p_h*(      K1 + p_h)
zalk_dic   = p_dictot * (znumer_dic/zdenom_dic)

! B(OH)3 - B(OH)4 : n=1, m=0
znumer_bor =       Kb
zdenom_bor =       Kb + p_h
zalk_bor   = p_bortot * (znumer_bor/zdenom_bor)

! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
znumer_po4 = 3._wp*K1p*K2p*K3p + p_h*(2._wp*K1p*K2p + p_h* K1p)
zdenom_po4 =       K1p*K2p*K3p + p_h*(      K1p*K2p + p_h*(K1p + p_h))
zalk_po4   = p_po4tot * (znumer_po4/zdenom_po4 - 1._wp) ! Zero level of H3PO4 = 1

! H4SiO4 - H3SiO4 : n=1, m=0
znumer_sil =       Ksi
zdenom_sil =       Ksi + p_h
zalk_sil   = p_siltot * (znumer_sil/zdenom_sil)

! NH4 - NH3 : n=1, m=0
!znumer_nh4 =       api1_nh4
!zdenom_nh4 =       api1_nh4 + p_h
!zalk_nh4   = p_nh4tot * (znumer_nh4/zdenom_nh4)
! Note: api1_nh4 = Knh4

! H2S - HS : n=1, m=0
!znumer_h2s =       api1_h2s
!zdenom_h2s =       api1_h2s + p_h
!zalk_h2s   = p_h2stot * (znumer_h2s/zdenom_h2s)
! Note: api1_h2s = Kh2s

! HSO4 - SO4 : n=1, m=1
znumer_so4 =       Ks
zdenom_so4 =       Ks + p_h
zalk_so4   = p_so4tot * (znumer_so4/zdenom_so4 - 1._wp)

! HF - F : n=1, m=1
znumer_flu =       Kf
zdenom_flu =       Kf + p_h
zalk_flu   = p_flutot * (znumer_flu/zdenom_flu - 1._wp)

! H2O - OH
zalk_wat   = Kw/p_h - p_h/aphscale

equation_at =    zalk_dic + zalk_bor + zalk_po4 + zalk_sil &
               + zalk_so4 + zalk_flu                       &
               + zalk_wat - p_alktot

IF(PRESENT(p_deriveqn)) THEN
   ! H2CO3 - HCO3 - CO3 : n=2
   zdnumer_dic = K1*K1*K2 + p_h*(4._wp*K1*K2               &
                          + p_h*       K1    )
   zdalk_dic   = -p_dictot*(zdnumer_dic/zdenom_dic**2)

   ! B(OH)3 - B(OH)4 : n=1
   zdnumer_bor = Kb
   zdalk_bor   = -p_bortot*(zdnumer_bor/zdenom_bor**2)

   ! H3PO4 - H2PO4 - HPO4 - PO4 : n=3
   zdnumer_po4 = K1p*K2p*K1p*K2p*K3p + p_h*(4._wp*K1p*K1p*K2p*K3p                &
                                     + p_h*(9._wp*K1p*K2p*K3p + K1p*K1p*K2p      &
                                     + p_h*(4._wp*K1p*K2p                        &
                                     + p_h*       K1p)))
   zdalk_po4   = -p_po4tot * (zdnumer_po4/zdenom_po4**2)

   ! H4SiO4 - H3SiO4 : n=1
   zdnumer_sil = Ksi
   zdalk_sil   = -p_siltot * (zdnumer_sil/zdenom_sil**2)

!  ! NH4 - NH3 : n=1
!  zdnumer_nh4 = Knh4
!  zdalk_nh4   = -p_nh4tot * (zdnumer_nh4/zdenom_nh4**2)

!  ! H2S - HS : n=1
!  zdnumer_h2s = api1_h2s
!  zdalk_h2s   = -p_h2stot * (zdnumer_h2s/zdenom_h2s**2)

   ! HSO4 - SO4 : n=1
   zdnumer_so4 = Ks
   zdalk_so4   = -p_so4tot * (zdnumer_so4/zdenom_so4**2)

   ! HF - F : n=1
   zdnumer_flu = Kf
   zdalk_flu   = -p_flutot * (zdnumer_flu/zdenom_flu**2)

!  p_deriveqn =   zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
!               + zdalk_nh4 + zdalk_h2s + zdalk_so4 + zdalk_flu &
!               - Kw/p_h**2 - 1._wp/aphscale
   p_deriveqn =   zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
                                        + zdalk_so4 + zdalk_flu &
                - Kw/p_h**2 - 1._wp/aphscale
ENDIF
RETURN
END FUNCTION equation_at

!===============================================================================

SUBROUTINE ahini_for_at(p_alkcb, p_dictot, p_bortot, K1, K2, Kb, p_hini)

! Subroutine returns the root for the 2nd order approximation of the
! DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
! around the local minimum, if it exists.

! Returns * 1E-03_wp if p_alkcb <= 0
!         * 1E-10_wp if p_alkcb >= 2*p_dictot + p_bortot
!         * 1E-07_wp if 0 < p_alkcb < 2*p_dictot + p_bortot
!                    and the 2nd order approximation does not have a solution

!USE MOD_CHEMCONST, ONLY : api1_dic, api2_dic, api1_bor

USE msingledouble
IMPLICIT NONE

! Argument variables
!--------------------
REAL(KIND=wp), INTENT(IN)   ::  p_alkcb, p_dictot, p_bortot
REAL(KIND=wp), INTENT(IN)   ::  K1, K2, Kb
REAL(KIND=wp), INTENT(OUT)  ::  p_hini

! Local variables
!-----------------
REAL(KIND=wp)  ::  zca, zba
REAL(KIND=wp)  ::  zd, zsqrtd, zhmin
REAL(KIND=wp)  ::  za2, za1, za0

IF (p_alkcb <= 0._wp) THEN
  p_hini = 1.e-3_wp
ELSEIF (p_alkcb >= (2._wp*p_dictot + p_bortot)) THEN
  p_hini = 1.e-10_wp
ELSE
  zca = p_dictot/p_alkcb
  zba = p_bortot/p_alkcb

  ! Coefficients of the cubic polynomial
  za2 = Kb*(1._wp - zba) + K1*(1._wp-zca)
  za1 = K1*Kb*(1._wp - zba - zca) + K1*K2*(1._wp - (zca+zca))
  za0 = K1*K2*Kb*(1._wp - zba - (zca+zca))
                                        ! Taylor expansion around the minimum
  zd = za2*za2 - 3._wp*za1              ! Discriminant of the quadratic equation
                                        ! for the minimum close to the root

  IF(zd > 0._wp) THEN                   ! If the discriminant is positive
    zsqrtd = SQRT(zd)
    IF(za2 < 0) THEN
      zhmin = (-za2 + zsqrtd)/3._wp
    ELSE
      zhmin = -za1/(za2 + zsqrtd)
    ENDIF
    p_hini = zhmin + SQRT(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)
  ELSE
    p_hini = 1.e-7_wp
  ENDIF

ENDIF
RETURN
END SUBROUTINE ahini_for_at

!===============================================================================

FUNCTION solve_at_general(p_alktot, p_dictot, p_bortot,                       &
                          p_po4tot, p_siltot,                                 &
                          p_so4tot, p_flutot,                                 &
                          K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,     &
                          p_hini,   p_val)

! Universal pH solver that converges from any given initial value,
! determines upper an lower bounds for the solution if required

USE msingledouble
IMPLICIT NONE
REAL(KIND=wp) :: SOLVE_AT_GENERAL

! Argument variables
!--------------------
REAL(KIND=wp), INTENT(IN)            :: p_alktot
REAL(KIND=wp), INTENT(IN)            :: p_dictot
REAL(KIND=wp), INTENT(IN)            :: p_bortot
REAL(KIND=wp), INTENT(IN)            :: p_po4tot
REAL(KIND=wp), INTENT(IN)            :: p_siltot
!REAL(KIND=wp), INTENT(IN)            :: p_nh4tot
!REAL(KIND=wp), INTENT(IN)            :: p_h2stot
REAL(KIND=wp), INTENT(IN)            :: p_so4tot
REAL(KIND=wp), INTENT(IN)            :: p_flutot
REAL(KIND=wp), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
REAL(KIND=wp), INTENT(IN)            :: K1p, K2p, K3p, Ksi
REAL(KIND=wp), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=wp), INTENT(OUT), OPTIONAL :: p_val

! Local variables
!-----------------
REAL(KIND=wp)  ::  zh_ini, zh, zh_prev, zh_lnfactor
REAL(KIND=wp)  ::  zalknw_inf, zalknw_sup
REAL(KIND=wp)  ::  zh_min, zh_max
REAL(KIND=wp)  ::  zdelta, zh_delta
REAL(KIND=wp)  ::  zeqn, zdeqndh, zeqn_absmin
REAL(KIND=wp)  ::  aphscale
LOGICAL        :: l_exitnow
REAL(KIND=wp), PARAMETER :: pz_exp_threshold = 1.0_wp

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
aphscale = 1._wp + p_so4tot/Ks

IF(PRESENT(p_hini)) THEN
   zh_ini = p_hini
ELSE
   CALL ahini_for_at(p_alktot, p_dictot, p_bortot, K1, K2, Kb, zh_ini)
ENDIF

 CALL anw_infsup(p_dictot, p_bortot,                                           &
                 p_po4tot, p_siltot,                                           &
                 p_so4tot, p_flutot,                                           &
                 zalknw_inf, zalknw_sup)

zdelta = (p_alktot-zalknw_inf)**2 + 4._wp*Kw/aphscale

IF(p_alktot >= zalknw_inf) THEN
   zh_min = 2._wp*Kw /( p_alktot-zalknw_inf + SQRT(zdelta) )
ELSE
   zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2._wp
ENDIF

zdelta = (p_alktot-zalknw_sup)**2 + 4._wp*Kw/aphscale

IF(p_alktot <= zalknw_sup) THEN
   zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2._wp
ELSE
   zh_max = 2._wp*Kw /( p_alktot-zalknw_sup + SQRT(zdelta) )
ENDIF

zh = MAX(MIN(zh_max, zh_ini), zh_min)
niter_atgen        = 0                 ! Reset counters of iterations
zeqn_absmin        = HUGE(1._wp)

DO
   IF(niter_atgen >= jp_maxniter_atgen) THEN
      zh = -1._wp
      EXIT
   ENDIF

   zh_prev = zh
   zeqn = equation_at(p_alktot, zh,       p_dictot, p_bortot,                  &
                      p_po4tot, p_siltot,                                      &
                      p_so4tot, p_flutot,                                      &
                      K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,          &
                      P_DERIVEQN = zdeqndh)

   ! Adapt bracketing interval
   IF(zeqn > 0._wp) THEN
      zh_min = zh_prev
   ELSEIF(zeqn < 0._wp) THEN
      zh_max = zh_prev
   ELSE
      ! zh is the root; unlikely but, one never knows
      EXIT
   ENDIF

   ! Now determine the next iterate zh
   niter_atgen = niter_atgen + 1

   IF(ABS(zeqn) >= 0.5_wp*zeqn_absmin) THEN
      ! if the function evaluation at the current point is
      ! not decreasing faster than with a bisection step (at least linearly)
      ! in absolute value take one bisection step on [ph_min, ph_max]
      ! ph_new = (ph_min + ph_max)/2d0
      !
      ! In terms of [H]_new:
      ! [H]_new = 10**(-ph_new)
      !         = 10**(-(ph_min + ph_max)/2d0)
      !         = SQRT(10**(-(ph_min + phmax)))
      !         = SQRT(zh_max * zh_min)
      zh = SQRT(zh_max * zh_min)
      zh_lnfactor = (zh - zh_prev)/zh_prev ! Required to test convergence below
   ELSE
      ! dzeqn/dpH = dzeqn/d[H] * d[H]/dpH
      !           = -zdeqndh * LOG(10) * [H]
      ! \Delta pH = -zeqn/(zdeqndh*d[H]/dpH) = zeqn/(zdeqndh*[H]*LOG(10))
      !
      ! pH_new = pH_old + \deltapH
      !
      ! [H]_new = 10**(-pH_new)
      !         = 10**(-pH_old - \Delta pH)
      !         = [H]_old * 10**(-zeqn/(zdeqndh*[H]_old*LOG(10)))
      !         = [H]_old * EXP(-LOG(10)*zeqn/(zdeqndh*[H]_old*LOG(10)))
      !         = [H]_old * EXP(-zeqn/(zdeqndh*[H]_old))

      zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

      IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
         zh          = zh_prev*EXP(zh_lnfactor)
      ELSE
         zh_delta    = zh_lnfactor*zh_prev
         zh          = zh_prev + zh_delta
      ENDIF

      IF( zh < zh_min ) THEN
         ! if [H]_new < [H]_min
         ! i.e., if ph_new > ph_max then
         ! take one bisection step on [ph_prev, ph_max]
         ! ph_new = (ph_prev + ph_max)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_max)/2d0)
         !         = SQRT(10**(-(ph_prev + phmax)))
         !         = SQRT([H]_old*10**(-ph_max))
         !         = SQRT([H]_old * zh_min)
         zh                = SQRT(zh_prev * zh_min)
         zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
      ENDIF

      IF( zh > zh_max ) THEN
         ! if [H]_new > [H]_max
         ! i.e., if ph_new < ph_min, then
         ! take one bisection step on [ph_min, ph_prev]
         ! ph_new = (ph_prev + ph_min)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_min)/2d0)
         !         = SQRT(10**(-(ph_prev + ph_min)))
         !         = SQRT([H]_old*10**(-ph_min))
         !         = SQRT([H]_old * zhmax)
         zh                = SQRT(zh_prev * zh_max)
         zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
      ENDIF
   ENDIF

   zeqn_absmin = MIN( ABS(zeqn), zeqn_absmin)

   ! Stop iterations once |\delta{[H]}/[H]| < rdel
   ! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
   ! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

   ! Alternatively:
   ! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
   !             ~ 1/LOG(10) * |\Delta [H]|/[H]
   !             < 1/LOG(10) * rdel

   ! Hence |zeqn/(zdeqndh*zh)| < rdel

   ! rdel <-- pp_rdel_ah_target

   l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

   IF(l_exitnow) EXIT
ENDDO

solve_at_general = zh

IF(PRESENT(p_val)) THEN
   IF(zh > 0._wp) THEN
      p_val = equation_at(p_alktot, zh,       p_dictot, p_bortot,              &
                          p_po4tot, p_siltot,                                  &
                          p_so4tot, p_flutot,                                  &
                          K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)
   ELSE
      p_val = HUGE(1._wp)
   ENDIF
ENDIF
RETURN
END FUNCTION solve_at_general

!===============================================================================

FUNCTION solve_at_general_sec(p_alktot, p_dictot, p_bortot,                   &
                              p_po4tot, p_siltot,                             &
                              p_so4tot, p_flutot,                             &
                              K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi, &
                              p_hini,   p_val)

! Universal pH solver that converges from any given initial value,
! determines upper an lower bounds for the solution if required

!USE MOD_CHEMCONST, ONLY: api1_wat, aphscale
USE msingledouble
IMPLICIT NONE
REAL(KIND=wp) :: SOLVE_AT_GENERAL_SEC

! Argument variables
REAL(KIND=wp), INTENT(IN)            :: p_alktot
REAL(KIND=wp), INTENT(IN)            :: p_dictot
REAL(KIND=wp), INTENT(IN)            :: p_bortot
REAL(KIND=wp), INTENT(IN)            :: p_po4tot
REAL(KIND=wp), INTENT(IN)            :: p_siltot
!REAL(KIND=wp), INTENT(IN)            :: p_nh4tot
!REAL(KIND=wp), INTENT(IN)            :: p_h2stot
REAL(KIND=wp), INTENT(IN)            :: p_so4tot
REAL(KIND=wp), INTENT(IN)            :: p_flutot
REAL(KIND=wp), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
REAL(KIND=wp), INTENT(IN)            :: K1p, K2p, K3p, Ksi
REAL(KIND=wp), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=wp), INTENT(OUT), OPTIONAL :: p_val

! Local variables
REAL(KIND=wp)  ::  zh_ini, zh, zh_1, zh_2, zh_delta
REAL(KIND=wp)  ::  zalknw_inf, zalknw_sup
REAL(KIND=wp)  ::  zh_min, zh_max
REAL(KIND=wp)  ::  zeqn, zeqn_1, zeqn_2, zeqn_absmin
REAL(KIND=wp)  ::  zdelta
REAL(KIND=wp)  ::  aphscale
LOGICAL        ::  l_exitnow

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
aphscale = 1._wp + p_so4tot/Ks

IF(PRESENT(p_hini)) THEN
   zh_ini = p_hini
ELSE
   CALL ahini_for_at(p_alktot, p_dictot, p_bortot, K1, K2, Kb, zh_ini)
ENDIF

 CALL anw_infsup(p_dictot, p_bortot,                                      &
                 p_po4tot, p_siltot,                                      &
                 p_so4tot, p_flutot,                                      &
                 zalknw_inf, zalknw_sup)

zdelta = (p_alktot-zalknw_inf)**2 + 4._wp*Kw/aphscale

IF(p_alktot >= zalknw_inf) THEN
   zh_min = 2._wp*Kw /( p_alktot-zalknw_inf + SQRT(zdelta) )
ELSE
   zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2._wp
ENDIF

zdelta = (p_alktot-zalknw_sup)**2 + 4._wp*Kw/aphscale

IF(p_alktot <= zalknw_sup) THEN
   zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2._wp
ELSE
   zh_max = 2._wp*Kw /( p_alktot-zalknw_sup + SQRT(zdelta) )
ENDIF

zh = MAX(MIN(zh_max, zh_ini), zh_min)
niter_atsec        = 0                 ! Reset counters of iterations

! Prepare the secant iterations: two initial (zh, zeqn) pairs are required
! We have the starting value, that needs to be completed by the evaluation
! of the equation value it produces.

! Complete the initial value with its equation evaluation
! (will take the role of the $n-2$ iterate at the first secant evaluation)

niter_atsec = 0                        ! zh_2 is the initial value;

zh_2   = zh
zeqn_2 = equation_at(p_alktot, zh_2,     p_dictot, p_bortot,                 &
                     p_po4tot, p_siltot,                                     &
                     p_so4tot, p_flutot,                                     &
                     K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)

zeqn_absmin        = ABS(zeqn_2)

! Adapt bracketing interval and heuristically set zh_1
IF(zeqn_2 < 0._wp) THEN
                                       ! If zeqn_2 < 0, then we adjust zh_max:
                                       ! we can be sure that zh_min < zh_2 < zh_max.
   zh_max = zh_2
                                       ! for zh_1, try 25% (0.1 pH units) below the current zh_max,
                                       ! but stay above SQRT(zh_min*zh_max), which would be equivalent
                                       ! to a bisection step on [pH@zh_min, pH@zh_max]
   zh_1   = MAX(zh_max/1.25_wp, SQRT(zh_min*zh_max))
ELSEIF(zeqn_2 > 0._wp) THEN
                                       ! If zeqn_2 < 0, then we adjust zh_min:
                                       ! we can be sure that zh_min < zh_2 < zh_max.
   zh_min = zh_2
                                       ! for zh_1, try 25% (0.1 pH units) above the current zh_min,
                                       ! but stay below SQRT(zh_min*zh_max) which would be equivalent
                                       ! to a bisection step on [pH@zh_min, pH@zh_max]
   zh_1   = MIN(zh_min*1.25_wp, SQRT(zh_min*zh_max))
ELSE ! we have got the root; unlikely, but one never knows
   solve_at_general_sec = zh_2
   IF(PRESENT(p_val)) p_val = zeqn_2
   RETURN
ENDIF

! We now have the first pair completed (zh_2, zeqn_2).
! Define the second one (zh_1, zeqn_1), which is also the first iterate.
! zh_1 has already been set above
niter_atsec = 1                        ! Update counter of iterations

zeqn_1 = equation_at(p_alktot, zh_1,       p_dictot, p_bortot,                 &
                     p_po4tot, p_siltot,                                       &
                     p_so4tot, p_flutot,                                       &
                     K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)

! Adapt bracketing interval: we know that zh_1 <= zh <= zh_max (if zeqn_1 > 0)
! or zh_min <= zh <= zh_1 (if zeqn_1 < 0), so this can always be done
IF(zeqn_1 > 0._wp) THEN
   zh_min = zh_1
ELSEIF(zeqn_1 < 0._wp) THEN
   zh_max = zh_1
ELSE ! zh_1 is the root
   solve_at_general_sec = zh_1
   IF(PRESENT(p_val)) p_val = zeqn_1
ENDIF

IF(ABS(zeqn_1) > zeqn_absmin) THEN     ! Swap zh_2 and zh_1 if ABS(zeqn_2) < ABS(zeqn_1)
                                       ! so that zh_2 and zh_1 lead to decreasing equation
                                       ! values (in absolute value)
   zh     = zh_1
   zeqn   = zeqn_1
   zh_1   = zh_2
   zeqn_1 = zeqn_2
   zh_2   = zh
   zeqn_2 = zeqn
ELSE
   zeqn_absmin = ABS(zeqn_1)
ENDIF

! Pre-calculate the first secant iterate (this is the second iterate)
niter_atsec = 2

zh_delta = -zeqn_1/((zeqn_2-zeqn_1)/(zh_2 - zh_1))
zh = zh_1 + zh_delta

! Make sure that zh_min < zh < zh_max (if not,
! bisect around zh_1 which is the best estimate)
IF (zh > zh_max) THEN                  ! this can only happen if zh_2 < zh_1
                                       ! and zeqn_2 > zeqn_1 > 0
   zh = SQRT(zh_1*zh_max)
ENDIF

IF (zh < zh_min) THEN                  ! this can only happen if zh_2 > zh_1
                                       ! and zeqn_2 < zeqn_1 < 0
   zh = SQRT(zh_1*zh_min)
ENDIF

DO
   IF(niter_atsec >= jp_maxniter_atsec) THEN
      zh = -1._wp
      EXIT
   ENDIF

   zeqn = equation_at(p_alktot, zh,       p_dictot, p_bortot,                  &
                      p_po4tot, p_siltot,                                      &
                      p_so4tot, p_flutot,                                      &
                      K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)

   ! Adapt bracketing interval: since initially, zh_min <= zh <= zh_max
   ! we are sure that zh will improve either bracket, depending on the sign
   ! of zeqn
   IF(zeqn > 0._wp) THEN
     zh_min = zh
   ELSEIF(zeqn < 0._wp) THEN
     zh_max = zh
   ELSE
     ! zh is the root
     EXIT
   ENDIF

   ! start calculation of next iterate
   niter_atsec = niter_atsec + 1

   zh_2   = zh_1
   zeqn_2 = zeqn_1
   zh_1   = zh
   zeqn_1 = zeqn

   IF(ABS(zeqn) >= 0.5_wp*zeqn_absmin) THEN
      ! if the function evaluation at the current point
      ! is not decreasing faster in absolute value than
      ! we may expect for a bisection step, then take
      ! one bisection step on [ph_min, ph_max]
      ! ph_new = (ph_min + ph_max)/2d0
      ! In terms of [H]_new:
      ! [H]_new = 10**(-ph_new)
      !         = 10**(-(ph_min + ph_max)/2d0)
      !         = SQRT(10**(-(ph_min + phmax)))
      !         = SQRT(zh_max * zh_min)
      zh                = SQRT(zh_max * zh_min)
      zh_delta           = zh - zh_1
   ELSE
      ! \Delta H = -zeqn_1*(h_2 - h_1)/(zeqn_2 - zeqn_1)
      ! H_new = H_1 + \Delta H
      zh_delta = -zeqn_1/((zeqn_2-zeqn_1)/(zh_2 - zh_1))
      zh       = zh_1 + zh_delta

      IF( zh < zh_min ) THEN
         ! if [H]_new < [H]_min
         ! i.e., if ph_new > ph_max then
         ! take one bisection step on [ph_prev, ph_max]
         ! ph_new = (ph_prev + ph_max)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_max)/2d0)
         !         = SQRT(10**(-(ph_prev + phmax)))
         !         = SQRT([H]_old*10**(-ph_max))
         !         = SQRT([H]_old * zh_min)
         zh                = SQRT(zh_1 * zh_min)
         zh_delta          = zh - zh_1
      ENDIF

      IF( zh > zh_max ) THEN
         ! if [H]_new > [H]_max
         ! i.e., if ph_new < ph_min, then
         ! take one bisection step on [ph_min, ph_prev]
         ! ph_new = (ph_prev + ph_min)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_min)/2d0)
         !         = SQRT(10**(-(ph_prev + ph_min)))
         !         = SQRT([H]_old*10**(-ph_min))
         !         = SQRT([H]_old * zhmax)
         zh                 = SQRT(zh_1 * zh_max)
         zh_delta           = zh - zh_1
      ENDIF
   ENDIF

   zeqn_absmin = MIN(ABS(zeqn), zeqn_absmin)

   ! Stop iterations once |([H]-[H_1])/[H_1]| < rdel
   l_exitnow = (ABS(zh_delta) < pp_rdel_ah_target*zh_1)

   IF(l_exitnow) EXIT
ENDDO

SOLVE_AT_GENERAL_SEC = zh

IF(PRESENT(p_val)) THEN
   IF(zh > 0._wp) THEN
      p_val = equation_at(p_alktot, zh,       p_dictot, p_bortot,              &
                          p_po4tot, p_siltot,                                  &
                          p_so4tot, p_flutot,                                  &
                          K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)
   ELSE
     p_val = HUGE(1._wp)
   ENDIF
ENDIF

RETURN
END FUNCTION SOLVE_AT_GENERAL_SEC

!===============================================================================

FUNCTION SOLVE_AT_FAST(p_alktot, p_dictot, p_bortot,                          &
                       p_po4tot, p_siltot,                                    &
                       p_so4tot, p_flutot,                                    &
                       K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,        &
                       p_hini,   p_val)

! Fast version of SOLVE_AT_GENERAL, without any bounds checking.

USE msingledouble
IMPLICIT NONE
REAL(KIND=wp) :: SOLVE_AT_FAST

! Argument variables
REAL(KIND=wp), INTENT(IN)            :: p_alktot
REAL(KIND=wp), INTENT(IN)            :: p_dictot
REAL(KIND=wp), INTENT(IN)            :: p_bortot
REAL(KIND=wp), INTENT(IN)            :: p_po4tot
REAL(KIND=wp), INTENT(IN)            :: p_siltot
!REAL(KIND=wp), INTENT(IN)            :: p_nh4tot
!REAL(KIND=wp), INTENT(IN)            :: p_h2stot
REAL(KIND=wp), INTENT(IN)            :: p_so4tot
REAL(KIND=wp), INTENT(IN)            :: p_flutot
REAL(KIND=wp), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
REAL(KIND=wp), INTENT(IN)            :: K1p, K2p, K3p, Ksi
REAL(KIND=wp), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=wp), INTENT(OUT), OPTIONAL :: p_val

! Local variables
REAL(KIND=wp)  ::  zh_ini, zh, zh_prev, zh_lnfactor
REAL(KIND=wp)  ::  zhdelta
REAL(KIND=wp)  ::  zeqn, zdeqndh
!REAL(KIND=wp)  ::  aphscale
LOGICAL        :: l_exitnow
REAL(KIND=wp), PARAMETER :: pz_exp_threshold = 1.0_wp

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
!aphscale = 1._wp + p_so4tot/Ks

IF(PRESENT(p_hini)) THEN
   zh_ini = p_hini
ELSE
   CALL AHINI_FOR_AT(p_alktot, p_dictot, p_bortot, K1, K2, Kb, zh_ini)
ENDIF
zh = zh_ini

niter_atfast    = 0                 ! Reset counters of iterations
DO
   niter_atfast = niter_atfast + 1
   IF(niter_atfast > jp_maxniter_atfast) THEN
      zh = -1._wp
      EXIT
   ENDIF

   zh_prev = zh

   zeqn = equation_at(p_alktot, zh,       p_dictot, p_bortot,                  &
                      p_po4tot, p_siltot,                                      &
                      p_so4tot, p_flutot,                                      &
                      K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,          &
                      P_DERIVEQN = zdeqndh)

   IF(zeqn == 0._wp) EXIT               ! zh is the root

   zh_lnfactor = -zeqn/(zdeqndh*zh_prev)
   IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
      zh          = zh_prev*EXP(zh_lnfactor)
   ELSE
      zhdelta     = zh_lnfactor*zh_prev
      zh          = zh_prev + zhdelta
   ENDIF

   ! Stop iterations once |\delta{[H]}/[H]| < rdel
   ! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
   ! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

   ! Alternatively:
   ! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
   !             ~ 1/LOG(10) * |\Delta [H]|/[H]
   !             < 1/LOG(10) * rdel

   ! Hence |zeqn/(zdeqndh*zh)| < rdel

   ! rdel <- pp_rdel_ah_target

   l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

   IF(l_exitnow) EXIT
ENDDO

SOLVE_AT_FAST = zh

IF(PRESENT(p_val)) THEN
   IF(zh > 0._wp) THEN
      p_val = equation_at(p_alktot, zh,       p_dictot, p_bortot,              &
                          p_po4tot, p_siltot,                                  &
                          p_so4tot, p_flutot,                                  &
                          K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)
   ELSE
      p_val = HUGE(1._wp)
   ENDIF
ENDIF

RETURN
END FUNCTION solve_at_fast
!===============================================================================

END MODULE mphsolvers


!> \file rho.f90
!! \BRIEF
!> Module with rho function - computes in situ density from S, T, P
MODULE mrho
CONTAINS
!> Function to compute in situ density from salinity (psu), in situ temperature (C), & pressure (bar)
FUNCTION rho(salt, temp, pbar)

  ! Compute in situ density from salinity (psu), in situ temperature (C), & pressure (bar)

  USE msingledouble
  IMPLICIT NONE

  !> salinity [psu]
  REAL(kind=r4) :: salt
  !> in situ temperature (C)
  REAL(kind=r4) :: temp
  !> pressure (bar) [Note units: this is NOT in decibars]
  REAL(kind=r4) :: pbar

  REAL(kind=r8) :: s, t, p
! REAL(kind=r8) :: t68
  REAL(kind=r8) :: X
  REAL(kind=r8) :: rhow, rho0
  REAL(kind=r8) :: a, b, c
  REAL(kind=r8) :: Ksbmw, Ksbm0, Ksbm
  REAL(kind=r8) :: drho

  REAL(kind=r4) :: rho

  !     Input arguments:
  !     -------------------------------------
  !     s  = salinity            [psu      (PSS-78) ]
  !     t  = in situ temperature [degree C (IPTS-68)]
  !     p  = pressure            [bar] !!!!  (not in [db]

  s = DBLE(salt)
  t = DBLE(temp)
  p = DBLE(pbar)

! Convert the temperature on today's "ITS 90" scale to older IPTS 68 scale
! (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
! According to Best-Practices guide, line above should be commented & 2 lines below should be uncommented
! Guide's answer of rho (s=35, t=25, p=0) = 1023.343 is for temperature given on ITPS-68 scale
! t68 = (T - 0.0002) / 0.99975
! X = t68
! Finally, don't do the ITS-90 to IPTS-68 conversion (T input var now already on IPTS-68 scale)
  X = T

! Density of pure water
  rhow = 999.842594d0 + 6.793952e-2_r8*X          &
       -9.095290e-3_r8*X*X + 1.001685e-4_r8*X**3  &
       -1.120083e-6_r8*X**4 + 6.536332e-9_r8*X**5

! Density of seawater at 1 atm, P=0
  A = 8.24493e-1_r8 - 4.0899e-3_r8*X                         &
       + 7.6438e-5_r8*X*X - 8.2467e-7_r8*X**3 + 5.3875e-9_r8*X**4
  B = -5.72466e-3_r8 + 1.0227e-4_r8*X - 1.6546e-6_r8*X*X
  C = 4.8314e-4_r8

  rho0 = rhow + A*S + B*S*SQRT(S) + C*S**2.0d0

! Secant bulk modulus of pure water
! The secant bulk modulus is the average change in pressure
! divided by the total change in volume per unit of initial volume.
  Ksbmw = 19652.21d0 + 148.4206d0*X - 2.327105d0*X*X &
       + 1.360477e-2_r8*X**3 - 5.155288e-5_r8*X**4

! Secant bulk modulus of seawater at 1 atm
  Ksbm0 = Ksbmw + S*( 54.6746d0 - 0.603459d0*X + 1.09987e-2_r8*X**2 &
       - 6.1670e-5_r8*X**3) &
       + S*SQRT(S)*( 7.944e-2_r8 + 1.6483e-2_r8*X - 5.3009e-4_r8*X**2)

! Secant bulk modulus of seawater at S,T,P
  Ksbm = Ksbm0 &
       + P*(3.239908d0 + 1.43713e-3_r8*X + 1.16092e-4_r8*X**2 - 5.77905e-7_r8*X**3) &
       + P*S*(2.2838e-3_r8 - 1.0981e-5_r8*X - 1.6078e-6_r8*X**2) &
       + P*S*SQRT(S)*1.91075e-4_r8 &
       + P*P*(8.50935e-5_r8 - 6.12293e-6_r8*X + 5.2787e-8_r8*X**2) &
       + P*P*S*(-9.9348e-7_r8 + 2.0816e-8_r8*X + 9.1697e-10_r8*X**2)

! Density of seawater at S,T,P
  drho = rho0/(1.0d0 - P/Ksbm)
  rho = REAL(drho)

  RETURN
END FUNCTION rho
END MODULE mrho


!> \file rhoinsitu.f90
!! \BRIEF
!> Module with rhoinsitu subroutine - compute in situ density from S, Tis, P
MODULE mrhoinsitu
CONTAINS
!>     Compute in situ density from salinity (psu), in situ temperature (C), & pressure (db).
!!     This subroutine is needed because rho is a function (using scalars not arrays)
SUBROUTINE rhoinsitu(salt, tempis, pdbar, N, rhois)

  !     Purpose:
  !     Compute in situ density from salinity (psu), in situ temperature (C), & pressure (db)
  !     Needed because rho is a function

  USE msingledouble
  USE mrho
  IMPLICIT NONE

  INTEGER :: N

! INPUT variables
  ! salt   = salinity [psu]
  ! tempis = in situ temperature [C]
  ! pdbar  = pressure [db]

  !> salinity [psu]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: salt
  !> in situ temperature [C]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: tempis
  !> pressure [db]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: pdbar
!f2py optional , depend(salt) :: n=len(salt)

! OUTPUT variables:
  ! rhois  = in situ density

  !> in situ density [kg/m3]
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: rhois

! Local variables
  INTEGER :: i

! REAL(kind=r4) ::  rho
! EXTERNAL rho

  DO i = 1,N
     rhois(i) = rho(salt(i), tempis(i), pdbar(i)/10.)
  END DO

  RETURN
END SUBROUTINE rhoinsitu
END MODULE mrhoinsitu


!> \file tis.f90
!! \BRIEF
!>    Module with tis subroutine - compute in situ T from S,T,P
MODULE mtis
CONTAINS
!>    Compute in situ temperature from arrays of potential temp, salinity, and pressure.
!!    This subroutine is needed because sw_temp is a function (using scalars not arrays)
SUBROUTINE tis(salt, tempot, press, pressref, N, tempis)
  !    Purpose:
  !    Compute in situ temperature from arrays of in situ temp, salinity, and pressure.
  !    Needed because sw_temp is a function

  USE msingledouble
  USE msw_temp
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> salinity [psu]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: salt
  !> potential temperature [C]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: tempot
  !> pressure [db]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: press
!f2py optional , depend(salt) :: n=len(salt)
  !> pressure reference level [db]
  REAL(kind=r4), INTENT(in) :: pressref

! OUTPUT variables:
  !> in situ temperature [C]
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: tempis

! REAL(kind=r8) :: dsalt, dtempis, dpress, dpressref
! REAL(kind=r8) :: dtempot

  INTEGER :: i

! REAL(kind=r8) :: sw_temp
! REAL(kind=r4) :: sw_temp
! EXTERNAL sw_temp

  DO i = 1,N
    !dsalt     = DBLE(salt(i))
    !dtempot   = DBLE(tempot(i))
    !dpress    = DBLE(press(i))
    !dpressref = DBLE(pressref)
    !dtempis   = sw_temp(dsalt, dtempot, dpress, dpressref)
    !tempis(i) = REAL(dtempis)

     tempis   = sw_temp(salt(i), tempot(i), press(i), pressref)
  END DO

  RETURN
END SUBROUTINE tis
END MODULE mtis


!> \file tpot.f90
!! \BRIEF
!>    Module with tpot subroutine - compute potential T from in situ T,S,P
MODULE mtpot
CONTAINS
!>    Compute potential temperature from arrays of in situ temp, salinity, and pressure.
!!    This subroutine is needed because sw_ptmp is a function (using scalars not arrays)
SUBROUTINE tpot(salt, tempis, press, pressref, N, tempot)
  !    Purpose:
  !    Compute potential temperature from arrays of in situ temp, salinity, and pressure.
  !    Needed because sw_ptmp is a function

  USE msingledouble
  USE msw_ptmp
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> salinity [psu]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: salt
  !> in situ temperature [C]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: tempis
  !> pressure [db]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: press
!f2py optional , depend(salt) :: n=len(salt)
  !> pressure reference level [db]
  REAL(kind=r4), INTENT(in) :: pressref

! OUTPUT variables:
  !> potential temperature [C] for pressref
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: tempot

  REAL(kind=r8) :: dsalt, dtempis, dpress, dpressref
  REAL(kind=r8) :: dtempot

  INTEGER :: i

! REAL(kind=r8) :: sw_ptmp
! EXTERNAL sw_ptmp

  DO i = 1,N
     dsalt     = DBLE(salt(i))
     dtempis   = DBLE(tempis(i))
     dpress    = DBLE(press(i))
     dpressref = DBLE(pressref)

     dtempot   = sw_ptmp(dsalt, dtempis, dpress, dpressref)

     tempot(i) = REAL(dtempot)
  END DO

  RETURN
END SUBROUTINE tpot
END MODULE mtpot



!> \file varsolver.f90
!! \BRIEF
!> Module with varsolver subroutine - solve for pH and other carbonate system variables
MODULE mvarsolver
CONTAINS
!>    Solve for pH and other carbonate system variables (with input from vars routine)
SUBROUTINE varsolver(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC,             &
                    temp, salt, ta, tc, pt, sit,                                 &
                    Bt, St, Ft,                                                  &
                    K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi,  &
                    Patm, Phydro_bar, rhodum, optGAS                             )

  !   Purpose: Solve for pH and other carbonate system variables (with input from vars routine)

  !     INPUT variables:
  !     ================
  !     temp    = in situ temperature [degrees C]
  !     ta      = total alkalinity                     in [eq/m^3] or [eq/kg]   based on optCON in calling routine (vars)
  !     tc      = dissolved inorganic carbon           in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     pt      = total dissolved inorganic phosphorus in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     sit     = total dissolved inorganic silicon    in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     Bt      = total dissolved inorganic boron      computed in calling routine (vars)
  !     St      = total dissolved inorganic sulfur     computed in calling routine (vars)
  !     Ft      = total dissolved inorganic fluorine   computed in calling routine (vars)
  !     K's     = K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi
  !     Patm    = atmospheric pressure [atm]
  !     Phydro_bar = hydrostatic pressure [bar]
  !     rhodum  = density factor as computed in calling routine  (vars)
  !     -----------
  !     optGAS: choose in situ vs. potential fCO2 and pCO2 (default optGAS = 'Pinsitu')
  !     ---------
  !       PRESSURE & T corrections for K0 and the fugacity coefficient (Cf)
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
  !     co2  = aqueous CO2 concentration in [mol/kg] or [mol/m^3] determined by rhodum (depends on optCON in calling routine)
  !     hco3 = bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] determined by rhodum
  !     co3  = carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] determined by rhodum
  !     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC = Omega for calcite, i.e., the   calcite saturation state

  USE msingledouble
  USE mphsolvers
  USE msw_ptmp

  IMPLICIT NONE

! Input variables
  !> <b>in situ temperature</b> [degrees C]
  REAL(kind=r8), INTENT(in) :: temp
  !> <b>salinity</b> [on the practical salinity scale, dimensionless]
  REAL(kind=r8), INTENT(in) :: salt
  !> total alkalinity in <b>[eq/m^3]</b> OR in <b>[eq/kg]</b>, depending on optCON in calling routine
  REAL(kind=r8), INTENT(in) :: ta
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  REAL(kind=r8), INTENT(in) :: tc
  !> phosphate concentration in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  REAL(kind=r8), INTENT(in) :: pt
  !> total dissolved inorganic silicon concentration in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  REAL(kind=r8), INTENT(in) :: sit
  !> total boron from either Uppstrom (1974) or Lee et al. (2010), depending on optB in calling routine
  REAL(kind=r8), INTENT(in) :: Bt
  !> total sulfate (Morris & Riley, 1966)
  REAL(kind=r8), INTENT(in) :: St
  !> total fluoride  (Riley, 1965)
  REAL(kind=r8), INTENT(in) :: Ft
  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  REAL(kind=r8), INTENT(in) :: K0
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), INTENT(in) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), INTENT(in) :: K2
  !> equilibrium constant for dissociation of boric acid
  REAL(kind=r8), INTENT(in) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  REAL(kind=r8), INTENT(in) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  REAL(kind=r8), INTENT(in) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride
  !! from Dickson and Riley (1979) or Perez and Fraga (1987), depending on optKf
  REAL(kind=r8), INTENT(in) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  REAL(kind=r8), INTENT(in) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  REAL(kind=r8), INTENT(in) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(in) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(in) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(in) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  REAL(kind=r8), INTENT(in) :: Ksi
  !> total atmospheric pressure <b>[atm]</b>
  REAL(kind=r8), INTENT(in) :: Patm
  !> total hydrostatic pressure <b>[bar]</b>
  REAL(kind=r8), INTENT(in) :: Phydro_bar
  !> density factor as computed incalling routine  (vars)
  REAL(kind=r8), INTENT(in) :: rhodum
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction)
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optGAS

! Output variables:
  !> pH on the <b>total scale</b>
  REAL(kind=r8), INTENT(out) :: ph
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(kind=r8), INTENT(out) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(kind=r8), INTENT(out) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=r8), INTENT(out) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=r8), INTENT(out) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=r8), INTENT(out) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=r8), INTENT(out) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(kind=r8), INTENT(out) :: OmegaC

! Local variables
  REAL(kind=r8) :: Phydro_atm, Ptot
  REAL(kind=r8) :: Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  REAL(kind=r8) :: tk, tk0
  real(kind=r8) :: temp68, tempot, tempot68
  REAL(kind=r8) :: Hinit, H
  REAL(kind=r8) :: HSO4, HF, HSI, HPO4
  REAL(kind=r8) :: ab, aw, ac
  REAL(kind=r8) :: cu, cb, cc
  REAL(kind=r8) :: Ca
! Array to pass optional arguments
  CHARACTER(7) :: opGAS

  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF

! Compute pH from constants and total concentrations
! - use SolveSAPHE v1.0.1 routines from Munhoven (2013, GMD) modified to use mocsy's Ks instead of its own
! 1) Compute best starting point for H+ calculation
  call ahini_for_at(ta, tc, Bt, K1, K2, Kb, Hinit)
! 2) Solve for H+ using above result as the initial H+ value
  H = solve_at_general(ta, tc, Bt,                                         &
                       pt,     sit,                                        &
                       St, Ft,                                             &
                       K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,     &
                       Hinit)
! 3) Calculate pH from from H+ concentration (mol/kg)
  IF (H > 0.d0) THEN
     pH = -1.*LOG10(H)
  ELSE
     pH = 1.e20_r8
  ENDIF

! Compute carbonate Alk (Ac) by difference: from total Alk and other Alk components
  HSO4 = St/(1.0d0 + Ks/(H/(1.0d0 + St/Ks)))
  HF = 1.0d0/(1.0d0 + Kf/H)
  HSI = 1.0d0/(1.0d0 + H/Ksi)
  HPO4 = (K1p*K2p*(H + 2.*K3p) - H**3) /                &
  (H**3 + K1p*H**2 + K1p*K2p*H + K1p*K2p*K3p)
  ab = Bt/(1.0d0 + H/Kb)
  aw = Kw/H - H/(1.0d0 + St/Ks)
  ac = ta + hso4 - sit*hsi - ab - aw + Ft*hf - pt*hpo4

! Calculate CO2*, HCO3-, & CO32- (in mol/kg soln) from Ct, Ac, H+, K1, & K2
  cu = (2.0d0 * tc - ac) / (2.0d0 + K1 / H)
  cb = K1 * cu / H
  cc = K2 * cb / H

! When optCON = 'mol/m3' in calling routine (vars), then:
! convert output var concentrations from mol/kg to mol/m^3
! e.g., for case when drho = 1028, multiply by [1.028 kg/L  x  1000 L/m^3])
  co2  = cu * rhodum
  hco3 = cb * rhodum
  co3  = cc * rhodum

! Determine CO2 fugacity [uatm]
! NOTE: equation just below requires CO2* in mol/kg
  fCO2 = cu * 1.e6_r8/K0

! Determine CO2 partial pressure from CO2 fugacity [uatm]
  tk = 273.15d0 + temp
  !Compute EITHER "potential pCO2" OR "in situ pCO2" (T and P used for calculations will differ)
  IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
     tk0 = tk                 !in situ temperature (K) for K0 calculation
     Ptot = Patm              !total pressure (in atm) = atmospheric pressure ONLY
  ELSEIF (trim(opGAS) == 'Ppot' .OR. trim(opGAS) == 'ppot') THEN
     !Use potential temperature and atmospheric pressure (water parcel adiabatically brought back to surface)
     !temp68 = (temp - 0.0002d0) / 0.99975d0          !temp = in situ T; temp68 is same converted to ITPS-68 scale
     !tempot68 = sw_ptmp(salt, temp68, Phydro_bar*10d0, 0.0d0) !potential temperature (C)
     !tempot   = 0.99975*tempot68 + 0.0002
     !tk0 = tempot + 273.15d0  !potential temperature (K) for fugacity coeff. calc as needed for potential fCO2 & pCO2
     tempot = sw_ptmp(salt, temp, Phydro_bar*10d0, 0.0d0) !potential temperature (C)
     tk0 = tempot + 273.15d0  !potential temperature (K) for fugacity coeff. calc as needed for potential fCO2 & pCO2
     Ptot = Patm              !total pressure (in atm) = atmospheric pressure ONLY
  ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
     !Use in situ temperature and total pressure
     tk0 = tk                             !in situ temperature (K) for fugacity coefficient calculation
     Phydro_atm = Phydro_bar / 1.01325d0  !convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
     Ptot = Patm + Phydro_atm            !total pressure (in atm) = atmospheric pressure + hydrostatic pressure
  ELSE
     PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
     STOP
  ENDIF

! Now that we have T and P in the right form, continue with calculation of fugacity coefficient (and pCO2)
  Rgas_atm = 82.05736_r8      ! (cm3 * atm) / (mol * K)  CODATA (2006)
! To compute fugcoeff, we need 3 other terms (B, Del, xc2) in addition to 3 others above (tk, Ptot, Rgas_atm)
  B = -1636.75d0 + 12.0408d0*tk0 - 0.0327957d0*(tk0*tk0) + 0.0000316528d0*(tk0*tk0*tk0)
  Del = 57.7d0 - 0.118d0*tk0
! "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
! x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
! Let's assume that xCO2 = fCO2. Resulting fugcoeff is identical to 8th digit after the decimal.
  xCO2approx = fCO2 * 1.e-6_r8
  IF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
!    xCO2approx = 400.0e-6_r8      !a simple test (gives about same result as seacarb for pCO2insitu)
!    approximate surface xCO2 ~ surface fCO2 (i.e., in situ fCO2 d by exponential pressure correction)
     xCO2approx = xCO2approx * exp( ((1-Ptot)*32.3_r8)/(82.05736_r8*tk0) )   ! of K0 press. correction, see Weiss (1974, equation 5)
  ENDIF
  xc2 = (1.0d0 - xCO2approx)**2
  fugcoeff = exp( Ptot*(B + 2.0d0*xc2*Del)/(Rgas_atm*tk0) )
  pCO2 = fCO2 / fugcoeff

! Determine Omega Calcite et Aragonite
! OmegaA = ((0.01028d0*salt/35.0d0)*cc) / Kspa
! OmegaC = ((0.01028d0*salt/35.0d0)*cc) / Kspc
! - see comments from Munhoven on the best value "0.02128" which differs slightly from the best practices guide (0.02127)
  Ca = (0.02128d0/40.078d0) * salt/1.80655d0
  OmegaA = (Ca*cc) / Kspa
  OmegaC = (Ca*cc) / Kspc

  RETURN
END SUBROUTINE varsolver
END MODULE mvarsolver


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


!> \file depth2press.f90
!! \BRIEF
!> Module with depth2press subroutine - converts depth to pressure
!! with Saunders (1981) formula
MODULE mdepth2press
CONTAINS
!>     Compute pressure [db] from depth [m] & latitude [degrees north].
!!     This subroutine is needed because p80 is a function (using scalars not arrays)
SUBROUTINE depth2press(depth, lat, pdbar, N)

  !     Purpose:
  !     Compute pressure [db] from depth [m] & latitude [degrees north].
  !     Needed because p80 is a function

  USE msingledouble
  USE mp80
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> depth [m]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: depth
  !> latitude [degrees]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: lat
!f2py optional , depend(depth) :: n=len(depth)

! OUTPUT variables:
  !> pressure [db]
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: pdbar

  !     Local variables
  INTEGER :: i

! REAL(kind=r4) ::  p80
! EXTERNAL p80

  DO i = 1,N
     pdbar(i) = p80(depth(i), lat(i))
  END DO

  RETURN
END SUBROUTINE depth2press
END MODULE mdepth2press


!> \file f2pCO2.f90
!! \BRIEF
!>    Module with f2pCO2 subroutine - compute pCO2 from fCO2, in situ T, atm pressure, hydrostatic pressure
MODULE mf2pCO2
CONTAINS
!>    Compute pCO2 from arrays of fCO2, in situ temp, atm pressure, & hydrostatic pressure
SUBROUTINE f2pCO2(fCO2, temp, Patm, p, N, pCO2)
  !    Purpose:
  !    Compute pCO2 from arrays of fCO2, in situ temp, atm pressure, & hydrostatic pressure

  USE msingledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> oceanic fugacity of CO2 [uatm]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: fCO2
  !> in situ temperature [C]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: temp
  !> atmospheric pressure [atm]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: Patm
  !> hydrostatic pressure [db]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: p
!f2py optional , depend(pCO2) :: n=len(pCO2)

! OUTPUT variables:
  !> oceanic partial pressure of CO2 [uatm]
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: pCO2

! LOCAL variables:
  REAL(kind=r8) :: dfCO2, dtemp, tk, dPatm, prb
  REAL(kind=r8) :: Ptot, Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  REAL(kind=r8) :: dpCO2

  INTEGER :: i

! REAL(kind=r8) :: sw_ptmp
! EXTERNAL sw_ptmp

  DO i = 1,N
     dfCO2     = DBLE(fCO2(i))
     dtemp     = DBLE(temp(i))
     dPatm     = DBLE(Patm(i))
     tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
     prb = DBLE(p(i)) / 10.0d0         !Pressure effect (prb is in bars)
     Ptot = dPatm + prb/1.01325d0      !Total pressure (atmospheric + hydrostatic) [atm]
     Rgas_atm = 82.05736_r8            !R in (cm3 * atm) / (mol * K)  from CODATA (2006)
!    To compute fugcoeff, we need 3 other terms (B, Del, xc2) as well as 3 others above (tk, Ptot, Rgas_atm)
     B = -1636.75d0 + 12.0408d0*tk - 0.0327957d0*(tk*tk) + 0.0000316528d0*(tk*tk*tk)
     Del = 57.7d0 - 0.118d0*tk
!    "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
!    x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
!    Let's assume that xCO2 = fCO2. Resulting fugcoeff is identical to 8th digit after the decimal.
     xCO2approx = dfCO2 * 1.e-6_r8
     xc2 = (1.0d0 - xCO2approx)**2
     fugcoeff = exp( Ptot*(B + 2.0d0*xc2*Del)/(Rgas_atm*tk) )
     dpCO2 = dfCO2 / fugcoeff
     pCO2(i) = REAL(dpCO2)
  END DO

  RETURN
END SUBROUTINE f2pCO2
END MODULE mf2pCO2


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
