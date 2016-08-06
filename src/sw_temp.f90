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
