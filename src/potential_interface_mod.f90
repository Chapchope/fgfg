module potential_interface_mod
  use types_mod
  implicit none
  private

  public :: potential_type
  public :: coulomb_2d_potential
  public :: create_coulomb_2d_potential

  !==============================================================================
  ! Abstract potential interface for easy replacement
  !==============================================================================
  type, abstract :: potential_type
   contains
     procedure(potential_value_interface), deferred :: value
     procedure(potential_force_interface), deferred :: force
     procedure(potential_moment_interface), deferred :: moment_coefficient
  end type potential_type

  abstract interface
     ! Compute potential at distance r
     pure function potential_value_interface(this, r, q) result(phi)
       import :: potential_type, dp
       class(potential_type), intent(in) :: this
       real(dp), intent(in) :: r, q
       real(dp) :: phi
     end function potential_value_interface

     ! Compute force at vector (rx, ry)
     pure subroutine potential_force_interface(this, rx, ry, q, fx, fy)
       import :: potential_type, dp
       class(potential_type), intent(in) :: this
       real(dp), intent(in) :: rx, ry, q
       real(dp), intent(out) :: fx, fy
     end subroutine potential_force_interface

     ! Multipole moment coefficient for order (n, m)
     ! Used in multipole expansion
     pure function potential_moment_interface(this, n, m, r) result(coeff)
       import :: potential_type, dp
       class(potential_type), intent(in) :: this
       integer, intent(in) :: n, m
       real(dp), intent(in) :: r
       real(dp) :: coeff
     end function potential_moment_interface
  end interface

  !==============================================================================
  ! 2D Coulomb potential: 1/r (z=0 plane)
  !==============================================================================
  type, extends(potential_type) :: coulomb_2d_potential
   contains
     procedure :: value => coulomb_2d_value
     procedure :: force => coulomb_2d_force
     procedure :: moment_coefficient => coulomb_2d_moment
  end type coulomb_2d_potential

contains

  !------------------------------------------------------------------------------
  ! Factory function to create Coulomb 2D potential
  !------------------------------------------------------------------------------
  function create_coulomb_2d_potential() result(pot)
    type(coulomb_2d_potential) :: pot
    ! Nothing to initialize for simple Coulomb (pot is default-initialized)
    ! This line prevents unused return value warning:
    pot = coulomb_2d_potential()
  end function create_coulomb_2d_potential

  !------------------------------------------------------------------------------
  ! Coulomb 2D potential implementations
  !------------------------------------------------------------------------------
  pure function coulomb_2d_value(this, r, q) result(phi)
    class(coulomb_2d_potential), intent(in) :: this
    real(dp), intent(in) :: r, q
    real(dp) :: phi

    if (r < TINY) then
      phi = 0.0_dp
    else
      phi = q / r
    end if
  end function coulomb_2d_value

  pure subroutine coulomb_2d_force(this, rx, ry, q, fx, fy)
    class(coulomb_2d_potential), intent(in) :: this
    real(dp), intent(in) :: rx, ry, q
    real(dp), intent(out) :: fx, fy

    real(dp) :: r2, r, r3

    r2 = rx*rx + ry*ry

    if (r2 < TINY*TINY) then
      fx = 0.0_dp
      fy = 0.0_dp
      return
    end if

    r = sqrt(r2)
    r3 = r * r2

    ! F = q * r / r^3
    fx = q * rx / r3
    fy = q * ry / r3
  end subroutine coulomb_2d_force

  pure function coulomb_2d_moment(this, n, m, r) result(coeff)
    class(coulomb_2d_potential), intent(in) :: this
    integer, intent(in) :: n, m
    real(dp), intent(in) :: r
    real(dp) :: coeff

    ! For Coulomb 1/r potential in 2D Cartesian expansion
    ! Multipole: r^n * (cos(m*theta) or sin(m*theta))
    ! Local: r^(-n-1) * (cos(m*theta) or sin(m*theta))

    if (r < TINY) then
      if (n == 0 .and. m == 0) then
        coeff = 1.0_dp
      else
        coeff = 0.0_dp
      end if
    else
      coeff = r**n
    end if
  end function coulomb_2d_moment

end module potential_interface_mod
