program test_grid_fmm
  !==============================================================================
  ! Test program for grid-based FMM with replaceable potential
  !
  ! Demonstrates:
  ! - Parametric expansion order
  ! - Accuracy validation against direct sum
  ! - Performance comparison
  ! - MD simulation example
  !==============================================================================
  use types_mod
  use potential_interface_mod
  use grid_fmm_mod
  use md_utils_mod
  implicit none

  integer :: N, nsteps, nframes
  real(dp) :: dt, Lbox
  real(dp), allocatable :: xs(:), ys(:), qs(:)
  real(dp), allocatable :: vxs(:), vys(:)
  real(dp), allocatable :: fx(:), fy(:), phi(:)
  real(dp), allocatable :: fx_dir(:), fy_dir(:), phi_dir(:)

  type(coulomb_2d_potential) :: potential
  type(grid_fmm_2d_params) :: fmm_params

  integer :: step, i, j
  real(dp) :: t0, t1, t_fmm, t_dir
  real(dp) :: rx, ry, r
  real(dp) :: max_err_phi, rms_err_phi, max_err_f, rms_err_f
  real(dp) :: err_phi, err_f
  integer :: imax_phi, imax_f
  integer :: dump_interval

  print *, "=============================================="
  print *, "  Grid-Based FMM for 2D Coulomb Potential"
  print *, "  with Parametric Order and Replaceable Pot."
  print *, "=============================================="
  print *, ""

  !============================================================================
  ! Configuration
  !============================================================================
  N = 10000             ! Number of particles
  dt = 0.001_dp         ! Time step
  Lbox = 1.0_dp         ! Box size
  nsteps = 5            ! MD steps
  nframes = 3           ! Output frames

  ! FMM parameters
  fmm_params%p_order = 8         ! Expansion order (higher = more accurate)
  fmm_params%ngrid_x = 40        ! Grid size
  fmm_params%ngrid_y = 40
  fmm_params%near_range = 2      ! Near-field range (cells)
  fmm_params%use_precompute = .true.
  fmm_params%box_margin = 0.1_dp

  print *, "Configuration:"
  print '(A,I0)', "  Number of particles (N): ", N
  print '(A,ES10.3)', "  Time step (dt): ", dt
  print '(A,F6.2)', "  Box size (Lbox): ", Lbox
  print '(A,I0)', "  MD steps: ", nsteps
  print '(A,I0)', "  Output frames: ", nframes
  print *, ""
  print *, "FMM Parameters:"
  print '(A,I0)', "  Expansion order (p): ", fmm_params%p_order
  print '(A,I0," x ",I0)', "  Grid size: ", fmm_params%ngrid_x, fmm_params%ngrid_y
  print '(A,I0)', "  Total grid cells: ", fmm_params%ngrid_x * fmm_params%ngrid_y
  print '(A,I0)', "  Near-field range: ", fmm_params%near_range, " cells"
  print *, ""

  !============================================================================
  ! Initialize particles
  !============================================================================
  allocate(xs(N), ys(N), qs(N))
  allocate(vxs(N), vys(N))
  allocate(fx(N), fy(N), phi(N))
  allocate(fx_dir(N), fy_dir(N), phi_dir(N))

  call random_seed()
  call random_number(xs)
  call random_number(ys)

  xs = (xs - 0.5_dp) * 2.0_dp * Lbox
  ys = (ys - 0.5_dp) * 2.0_dp * Lbox
  qs = 1.0_dp  ! All charges equal to 1
  vxs = 0.0_dp
  vys = 0.0_dp

  ! Create potential
  potential = create_coulomb_2d_potential()
  print *, "Potential: Coulomb 2D (1/r)"
  print *, ""

  !============================================================================
  ! Initial validation: FMM vs Direct Sum
  !============================================================================
  print *, "Validation Test (Step 0):"
  print *, "----------------------------------------"

  ! FMM computation
  call cpu_time(t0)
  call grid_fmm_2d_compute(N, xs, ys, qs, potential, fmm_params, phi, fx, fy)
  call cpu_time(t1)
  t_fmm = t1 - t0

  print '(A,F10.6,A)', "  FMM time: ", t_fmm, " s"

  ! Direct sum (for validation)
  print *, "  Computing direct sum..."
  call cpu_time(t0)
  call direct_sum_2d(N, xs, ys, qs, potential, phi_dir, fx_dir, fy_dir)
  call cpu_time(t1)
  t_dir = t1 - t0

  print '(A,F10.6,A)', "  Direct sum time: ", t_dir, " s"
  print '(A,F8.2,A)', "  Speedup: ", t_dir / t_fmm, "x"
  print *, ""

  !============================================================================
  ! Error analysis
  !============================================================================
  print *, "Accuracy Analysis:"
  print *, "----------------------------------------"

  ! Potential errors
  max_err_phi = 0.0_dp
  rms_err_phi = 0.0_dp
  imax_phi = 1

  do i = 1, N
    err_phi = abs(phi(i) - phi_dir(i))
    rms_err_phi = rms_err_phi + err_phi**2

    if (err_phi > max_err_phi) then
      max_err_phi = err_phi
      imax_phi = i
    end if
  end do
  rms_err_phi = sqrt(rms_err_phi / real(N, dp))

  ! Force errors
  max_err_f = 0.0_dp
  rms_err_f = 0.0_dp
  imax_f = 1

  do i = 1, N
    err_f = sqrt((fx(i) - fx_dir(i))**2 + (fy(i) - fy_dir(i))**2)
    rms_err_f = rms_err_f + err_f**2

    if (err_f > max_err_f) then
      max_err_f = err_f
      imax_f = i
    end if
  end do
  rms_err_f = sqrt(rms_err_f / real(N, dp))

  print *, "Potential:"
  print '(A,ES12.4)', "  Max absolute error: ", max_err_phi
  print '(A,I0)', "    (at particle ", imax_phi, ")"
  print '(A,ES12.4)', "  RMS error: ", rms_err_phi
  print *, ""

  print *, "Force:"
  print '(A,ES12.4)', "  Max absolute error: ", max_err_f
  print '(A,I0)', "    (at particle ", imax_f, ")"
  print '(A,ES12.4)', "  RMS error: ", rms_err_f
  print *, ""

  ! Sample comparison (first 5 particles)
  print *, "Sample values (first 5 particles):"
  print *, "  i  |    phi (FMM)    |   phi (Direct)  |     fx (FMM)    |   fx (Direct)"
  print *, "-----+-----------------+-----------------+-----------------+----------------"
  do i = 1, min(5, N)
    print '(I5," | ",ES14.6," | ",ES14.6," | ",ES14.6," | ",ES14.6)', &
          i, phi(i), phi_dir(i), fx(i), fx_dir(i)
  end do
  print *, ""

  !============================================================================
  ! MD Simulation
  !============================================================================
  print *, "=============================================="
  print *, "Starting MD simulation..."
  print *, "=============================================="
  print *, ""

  dump_interval = max(1, nsteps / (nframes - 1))

  call write_frame_csv(0, N, xs, ys)
  print '(A)', "  Wrote frame_00000.csv"

  do step = 1, nsteps
    print '(A,I0,A,I0)', "Step ", step, " / ", nsteps

    ! Compute forces using FMM
    call cpu_time(t0)
    call grid_fmm_2d_compute(N, xs, ys, qs, potential, fmm_params, phi, fx, fy)
    call cpu_time(t1)
    print '(A,F10.6,A)', "  FMM time: ", t1 - t0, " s"

    ! Integrate
    call md_step_euler(N, xs, ys, vxs, vys, fx, fy, dt)

    ! Apply boundary conditions
    call apply_reflect_bc(N, xs, ys, vxs, vys, Lbox)

    ! Write output
    if (mod(step, dump_interval) == 0 .or. step == nsteps) then
      call write_frame_csv(step, N, xs, ys)
      print '(A,I5.5,A)', "  Wrote frame_", step, ".csv"
    end if

    print *, ""
  end do

  print *, "=============================================="
  print *, "Simulation completed successfully!"
  print *, "=============================================="

  ! Cleanup
  deallocate(xs, ys, qs, vxs, vys)
  deallocate(fx, fy, phi)
  deallocate(fx_dir, fy_dir, phi_dir)

contains

  !----------------------------------------------------------------------------
  ! Direct sum for validation
  !----------------------------------------------------------------------------
  subroutine direct_sum_2d(N, xs, ys, qs, potential, phi, fx, fy)
    integer, intent(in) :: N
    real(dp), intent(in) :: xs(N), ys(N), qs(N)
    class(potential_type), intent(in) :: potential
    real(dp), intent(out) :: phi(N), fx(N), fy(N)

    integer :: i, j
    real(dp) :: rx, ry, r, phi_ij, fx_ij, fy_ij

    phi = 0.0_dp
    fx = 0.0_dp
    fy = 0.0_dp

    do i = 1, N
      do j = 1, N
        if (i == j) cycle

        rx = xs(i) - xs(j)
        ry = ys(i) - ys(j)
        r = sqrt(rx*rx + ry*ry)

        if (r < TINY) cycle

        phi_ij = potential%value(r, qs(j))
        call potential%force(rx, ry, qs(j), fx_ij, fy_ij)

        phi(i) = phi(i) + phi_ij
        fx(i) = fx(i) + fx_ij
        fy(i) = fy(i) + fy_ij
      end do
    end do

  end subroutine direct_sum_2d

end program test_grid_fmm
