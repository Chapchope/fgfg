program demo_quick
  !==============================================================================
  ! Quick demo of grid-based FMM
  !==============================================================================
  use types_mod
  use potential_interface_mod
  use grid_fmm_mod
  implicit none

  integer :: N, i
  real(dp), allocatable :: xs(:), ys(:), qs(:)
  real(dp), allocatable :: phi(:), fx(:), fy(:)
  type(coulomb_2d_potential) :: potential
  type(grid_fmm_2d_params) :: fmm_params

  real(dp) :: t0, t1

  print *, "=============================================="
  print *, "  Quick Demo: Grid-Based FMM"
  print *, "=============================================="
  print *, ""

  ! Configuration
  N = 1000
  fmm_params%p_order = 6
  fmm_params%ngrid_x = 20
  fmm_params%ngrid_y = 20
  fmm_params%near_range = 2
  fmm_params%use_precompute = .true.

  print *, "Configuration:"
  print '(A,I0)', "  Particles: ", N
  print '(A,I0)', "  Expansion order: ", fmm_params%p_order
  print '(A,I0," x ",I0)', "  Grid: ", fmm_params%ngrid_x, fmm_params%ngrid_y
  print *, ""

  ! Initialize
  allocate(xs(N), ys(N), qs(N))
  allocate(phi(N), fx(N), fy(N))

  call random_seed()
  call random_number(xs)
  call random_number(ys)

  xs = (xs - 0.5_dp) * 2.0_dp
  ys = (ys - 0.5_dp) * 2.0_dp
  qs = 1.0_dp

  potential = create_coulomb_2d_potential()

  ! Run FMM
  print *, "Running FMM..."
  call cpu_time(t0)
  call grid_fmm_2d_compute(N, xs, ys, qs, potential, fmm_params, phi, fx, fy)
  call cpu_time(t1)

  print '(A,F10.6,A)', "FMM time: ", t1 - t0, " s"
  print *, ""

  ! Show sample results
  print *, "Sample results (first 5 particles):"
  print *, "  i  |      phi        |       fx        |       fy"
  print *, "-----+-----------------+-----------------+----------------"
  do i = 1, min(5, N)
    print '(I5," | ",ES14.6," | ",ES14.6," | ",ES14.6)', i, phi(i), fx(i), fy(i)
  end do
  print *, ""

  print *, "Demo completed successfully!"
  print *, "=============================================="

  deallocate(xs, ys, qs, phi, fx, fy)

end program demo_quick
