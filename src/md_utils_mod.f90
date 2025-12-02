module md_utils_mod
  !==============================================================================
  ! Molecular dynamics utilities
  !==============================================================================
  use types_mod
  implicit none
  private

  public :: md_step_euler
  public :: apply_periodic_bc
  public :: apply_reflect_bc
  public :: write_frame_csv
  public :: compute_kinetic_energy
  public :: compute_total_energy

contains

  !----------------------------------------------------------------------------
  ! Euler integration step
  !----------------------------------------------------------------------------
  subroutine md_step_euler(N, xs, ys, vxs, vys, fx, fy, dt)
    integer, intent(in) :: N
    real(dp), intent(inout) :: xs(N), ys(N), vxs(N), vys(N)
    real(dp), intent(in) :: fx(N), fy(N)
    real(dp), intent(in) :: dt

    integer :: i

    do i = 1, N
      vxs(i) = vxs(i) + fx(i) * dt
      vys(i) = vys(i) + fy(i) * dt
      xs(i) = xs(i) + vxs(i) * dt
      ys(i) = ys(i) + vys(i) * dt
    end do

  end subroutine md_step_euler

  !----------------------------------------------------------------------------
  ! Apply periodic boundary conditions
  !----------------------------------------------------------------------------
  subroutine apply_periodic_bc(N, xs, ys, Lx, Ly, xmin, ymin)
    integer, intent(in) :: N
    real(dp), intent(inout) :: xs(N), ys(N)
    real(dp), intent(in) :: Lx, Ly, xmin, ymin

    integer :: i
    real(dp) :: xmax, ymax

    xmax = xmin + Lx
    ymax = ymin + Ly

    do i = 1, N
      ! Wrap in x
      do while (xs(i) < xmin)
        xs(i) = xs(i) + Lx
      end do
      do while (xs(i) >= xmax)
        xs(i) = xs(i) - Lx
      end do

      ! Wrap in y
      do while (ys(i) < ymin)
        ys(i) = ys(i) + Ly
      end do
      do while (ys(i) >= ymax)
        ys(i) = ys(i) - Ly
      end do
    end do

  end subroutine apply_periodic_bc

  !----------------------------------------------------------------------------
  ! Apply reflecting boundary conditions
  !----------------------------------------------------------------------------
  subroutine apply_reflect_bc(N, xs, ys, vxs, vys, L)
    integer, intent(in) :: N
    real(dp), intent(inout) :: xs(N), ys(N), vxs(N), vys(N)
    real(dp), intent(in) :: L

    integer :: i

    do i = 1, N
      if (xs(i) > L) then
        xs(i) = 2.0_dp * L - xs(i)
        vxs(i) = -vxs(i)
      else if (xs(i) < -L) then
        xs(i) = -2.0_dp * L - xs(i)
        vxs(i) = -vxs(i)
      end if

      if (ys(i) > L) then
        ys(i) = 2.0_dp * L - ys(i)
        vys(i) = -vys(i)
      else if (ys(i) < -L) then
        ys(i) = -2.0_dp * L - ys(i)
        vys(i) = -vys(i)
      end if
    end do

  end subroutine apply_reflect_bc

  !----------------------------------------------------------------------------
  ! Write frame to CSV file
  !----------------------------------------------------------------------------
  subroutine write_frame_csv(step, N, xs, ys, filename_prefix)
    integer, intent(in) :: step, N
    real(dp), intent(in) :: xs(N), ys(N)
    character(len=*), intent(in), optional :: filename_prefix

    character(len=256) :: fname
    character(len=64) :: prefix
    integer :: i, unit

    if (present(filename_prefix)) then
      prefix = trim(filename_prefix)
    else
      prefix = 'frame'
    end if

    write(fname, '(A,A,I5.5,A)') trim(prefix), '_', step, '.csv'

    open(newunit=unit, file=fname, status='replace', action='write')
    write(unit, '(A)') 'x,y'
    do i = 1, N
      write(unit, '(ES16.8,",",ES16.8)') xs(i), ys(i)
    end do
    close(unit)

  end subroutine write_frame_csv

  !----------------------------------------------------------------------------
  ! Compute kinetic energy
  !----------------------------------------------------------------------------
  function compute_kinetic_energy(N, vxs, vys, masses) result(KE)
    integer, intent(in) :: N
    real(dp), intent(in) :: vxs(N), vys(N)
    real(dp), intent(in), optional :: masses(N)
    real(dp) :: KE

    integer :: i
    real(dp) :: m

    KE = 0.0_dp

    if (present(masses)) then
      do i = 1, N
        KE = KE + 0.5_dp * masses(i) * (vxs(i)**2 + vys(i)**2)
      end do
    else
      do i = 1, N
        KE = KE + 0.5_dp * (vxs(i)**2 + vys(i)**2)
      end do
    end if

  end function compute_kinetic_energy

  !----------------------------------------------------------------------------
  ! Compute total energy (kinetic + potential)
  !----------------------------------------------------------------------------
  function compute_total_energy(N, vxs, vys, phi, masses) result(E_total)
    integer, intent(in) :: N
    real(dp), intent(in) :: vxs(N), vys(N), phi(N)
    real(dp), intent(in), optional :: masses(N)
    real(dp) :: E_total

    real(dp) :: KE, PE
    integer :: i

    KE = compute_kinetic_energy(N, vxs, vys, masses)

    ! Potential energy (factor of 0.5 to avoid double counting)
    PE = 0.0_dp
    do i = 1, N
      PE = PE + 0.5_dp * phi(i)
    end do

    E_total = KE + PE

  end function compute_total_energy

end module md_utils_mod
