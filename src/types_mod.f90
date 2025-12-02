module types_mod
  implicit none

  ! Double precision kind
  integer, parameter :: dp = kind(1.0d0)

  ! Mathematical constants
  real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: SQRT_PI = 1.772453850905516027298167483341145182798_dp
  real(dp), parameter :: TWO_PI = 6.283185307179586476925286766559005768394_dp

  ! Numerical constants
  real(dp), parameter :: TINY = 1.0e-14_dp
  real(dp), parameter :: SMALL = 1.0e-12_dp

end module types_mod
