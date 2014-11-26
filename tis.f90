!>    Module with tis subroutine - compute in situ T from S,T,P
MODULE mtis
CONTAINS
!>    Compute in situ temperature from arrays of potential temp, salinity, and pressure.
!!    This subroutine is needed because sw_temp is a function (cannot accept arrays)
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
