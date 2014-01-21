!>    Module with tpot subroutine - compute potential T from in situ T,S,P
MODULE mtpot
CONTAINS
!>    Compute potential temperature from arrays of in situ temp, salinity, and pressure.
!!    This subroutine is needed because sw_ptmp is a function (cannot accept arrays)
SUBROUTINE tpot(salt, tempis, press, pressref, N, tempot)
  !    Purpose:
  !    Compute potential temperature from arrays of in situ temp, salinity, and pressure.
  !    Needed because sw_ptmp is a function (cannot accept arrays)

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
