module grid_fmm_mod
  !==============================================================================
  ! Grid-based Fast Multipole Method for 2D problems
  !
  ! Features:
  ! - Uniform grid structure (no tree building overhead)
  ! - Parametric expansion order (user-selectable accuracy)
  ! - Replaceable potential interface
  ! - O(N) complexity for large N
  ! - High accuracy (10^-10 achievable with high p)
  !
  ! Algorithm:
  ! 1. P2M: Particles to multipole moments (bottom-up)
  ! 2. M2M: Multipole to multipole (coarse grid, optional)
  ! 3. M2L: Multipole to local (far-field interactions)
  ! 4. L2L: Local to local (top-down, optional)
  ! 5. L2P + P2P: Local evaluation + near-field direct
  !==============================================================================
  use types_mod
  use potential_interface_mod
  use cartesian_multipole_mod
  implicit none
  private

  public :: grid_fmm_2d_compute
  public :: grid_fmm_2d_params

  !------------------------------------------------------------------------------
  ! Grid FMM parameters
  !------------------------------------------------------------------------------
  type :: grid_fmm_2d_params
     integer :: p_order = 6           ! Multipole expansion order (default: 6)
     integer :: ngrid_x = 32          ! Grid cells in x direction
     integer :: ngrid_y = 32          ! Grid cells in y direction
     integer :: near_range = 2        ! Near-field range (cells)
     logical :: use_precompute = .true.  ! Precompute factorials
     real(dp) :: box_margin = 0.1_dp  ! Margin around particles (fraction)
  end type grid_fmm_2d_params

  !------------------------------------------------------------------------------
  ! Grid cell structure
  !------------------------------------------------------------------------------
  type :: grid_cell_2d
     real(dp) :: xc, yc                       ! Cell center
     integer :: np                            ! Number of particles
     integer, allocatable :: plist(:)         ! Particle indices
     real(dp), allocatable :: M(:,:)          ! Multipole moments (0:p, 0:p)
     real(dp), allocatable :: L(:,:)          ! Local expansion (0:p, 0:p)
  end type grid_cell_2d

contains

  !============================================================================
  ! Main FMM computation routine
  !============================================================================
  subroutine grid_fmm_2d_compute(N, xs, ys, qs, potential, params, phi, fx, fy)
    !--------------------------------------------------------------------------
    ! Input:
    !   N         - Number of particles
    !   xs, ys    - Particle coordinates (z=0 assumed)
    !   qs        - Particle charges
    !   potential - Potential object (must extend potential_type)
    !   params    - FMM parameters
    !
    ! Output:
    !   phi       - Potential at each particle
    !   fx, fy    - Force at each particle
    !--------------------------------------------------------------------------
    integer, intent(in) :: N
    real(dp), intent(in) :: xs(N), ys(N), qs(N)
    class(potential_type), intent(in) :: potential
    type(grid_fmm_2d_params), intent(in) :: params
    real(dp), intent(out) :: phi(N), fx(N), fy(N)

    type(grid_cell_2d), allocatable :: grid(:,:)
    real(dp) :: xmin, xmax, ymin, ymax, Lx, Ly, hx, hy
    integer :: ngx, ngy, p_order
    integer :: ix, iy

    ! Extract parameters
    ngx = params%ngrid_x
    ngy = params%ngrid_y
    p_order = params%p_order

    ! Precompute factorials if requested
    if (params%use_precompute) then
      call precompute_factorials(2 * p_order + 2)
    end if

    ! Determine domain bounds
    call compute_domain_bounds(N, xs, ys, params%box_margin, &
                              xmin, xmax, ymin, ymax, Lx, Ly)

    hx = Lx / real(ngx, dp)
    hy = Ly / real(ngy, dp)

    ! Allocate and initialize grid
    allocate(grid(ngx, ngy))
    call initialize_grid(grid, ngx, ngy, xmin, ymin, hx, hy, p_order)

    ! Assign particles to cells
    call assign_particles_to_grid(N, xs, ys, grid, ngx, ngy, xmin, ymin, hx, hy)

    ! Step 1: P2M - Compute multipole expansions
    call p2m_step(grid, ngx, ngy, xs, ys, qs, p_order)

    ! Step 2: M2L - Far-field interactions (multipole to local)
    call m2l_step(grid, ngx, ngy, p_order, params%near_range)

    ! Step 3: L2P + P2P - Evaluate local + near-field direct
    phi = 0.0_dp
    fx = 0.0_dp
    fy = 0.0_dp

    call evaluate_step(grid, ngx, ngy, N, xs, ys, qs, potential, &
                      p_order, params%near_range, phi, fx, fy)

    ! Cleanup
    call cleanup_grid(grid, ngx, ngy)
    deallocate(grid)

  end subroutine grid_fmm_2d_compute

  !============================================================================
  ! Helper subroutines
  !============================================================================

  !----------------------------------------------------------------------------
  ! Compute domain bounds with margin
  !----------------------------------------------------------------------------
  subroutine compute_domain_bounds(N, xs, ys, margin, xmin, xmax, ymin, ymax, Lx, Ly)
    integer, intent(in) :: N
    real(dp), intent(in) :: xs(N), ys(N), margin
    real(dp), intent(out) :: xmin, xmax, ymin, ymax, Lx, Ly

    real(dp) :: xcen, ycen

    xmin = minval(xs)
    xmax = maxval(xs)
    ymin = minval(ys)
    ymax = maxval(ys)

    Lx = xmax - xmin
    Ly = ymax - ymin

    ! Add margin
    Lx = Lx * (1.0_dp + margin)
    Ly = Ly * (1.0_dp + margin)

    ! Center domain
    xcen = 0.5_dp * (xmin + xmax)
    ycen = 0.5_dp * (ymin + ymax)

    xmin = xcen - 0.5_dp * Lx
    xmax = xcen + 0.5_dp * Lx
    ymin = ycen - 0.5_dp * Ly
    ymax = ycen + 0.5_dp * Ly

  end subroutine compute_domain_bounds

  !----------------------------------------------------------------------------
  ! Initialize grid structure
  !----------------------------------------------------------------------------
  subroutine initialize_grid(grid, ngx, ngy, xmin, ymin, hx, hy, p_order)
    type(grid_cell_2d), intent(inout) :: grid(:,:)
    integer, intent(in) :: ngx, ngy, p_order
    real(dp), intent(in) :: xmin, ymin, hx, hy

    integer :: ix, iy

    do iy = 1, ngy
      do ix = 1, ngx
        grid(ix, iy)%xc = xmin + (real(ix, dp) - 0.5_dp) * hx
        grid(ix, iy)%yc = ymin + (real(iy, dp) - 0.5_dp) * hy
        grid(ix, iy)%np = 0

        allocate(grid(ix, iy)%M(0:p_order, 0:p_order))
        allocate(grid(ix, iy)%L(0:p_order, 0:p_order))

        grid(ix, iy)%M = 0.0_dp
        grid(ix, iy)%L = 0.0_dp
      end do
    end do

  end subroutine initialize_grid

  !----------------------------------------------------------------------------
  ! Assign particles to grid cells
  !----------------------------------------------------------------------------
  subroutine assign_particles_to_grid(N, xs, ys, grid, ngx, ngy, xmin, ymin, hx, hy)
    integer, intent(in) :: N, ngx, ngy
    real(dp), intent(in) :: xs(N), ys(N)
    type(grid_cell_2d), intent(inout) :: grid(:,:)
    real(dp), intent(in) :: xmin, ymin, hx, hy

    integer :: ip, ix, iy, i
    integer, allocatable :: counts(:,:)
    integer, allocatable :: temp_lists(:,:,:)

    ! Count particles per cell
    allocate(counts(ngx, ngy))
    counts = 0

    do ip = 1, N
      ix = int((xs(ip) - xmin) / hx) + 1
      iy = int((ys(ip) - ymin) / hy) + 1

      ! Clamp to bounds
      ix = max(1, min(ngx, ix))
      iy = max(1, min(ngy, iy))

      counts(ix, iy) = counts(ix, iy) + 1
    end do

    ! Allocate temporary storage
    allocate(temp_lists(ngx, ngy, maxval(counts)))

    ! Reset counts for indexing
    do iy = 1, ngy
      do ix = 1, ngx
        grid(ix, iy)%np = counts(ix, iy)
        if (counts(ix, iy) > 0) then
          allocate(grid(ix, iy)%plist(counts(ix, iy)))
        else
          allocate(grid(ix, iy)%plist(0))
        end if
      end do
    end do
    counts = 0

    ! Assign particles
    do ip = 1, N
      ix = int((xs(ip) - xmin) / hx) + 1
      iy = int((ys(ip) - ymin) / hy) + 1

      ix = max(1, min(ngx, ix))
      iy = max(1, min(ngy, iy))

      counts(ix, iy) = counts(ix, iy) + 1
      temp_lists(ix, iy, counts(ix, iy)) = ip
    end do

    ! Copy to grid
    do iy = 1, ngy
      do ix = 1, ngx
        if (grid(ix, iy)%np > 0) then
          grid(ix, iy)%plist = temp_lists(ix, iy, 1:grid(ix, iy)%np)
        end if
      end do
    end do

    deallocate(counts, temp_lists)

  end subroutine assign_particles_to_grid

  !----------------------------------------------------------------------------
  ! P2M step: Compute multipole moments for each cell
  !----------------------------------------------------------------------------
  subroutine p2m_step(grid, ngx, ngy, xs, ys, qs, p_order)
    type(grid_cell_2d), intent(inout) :: grid(:,:)
    integer, intent(in) :: ngx, ngy, p_order
    real(dp), intent(in) :: xs(:), ys(:), qs(:)

    integer :: ix, iy, np
    real(dp), allocatable :: xs_cell(:), ys_cell(:), qs_cell(:)

    do iy = 1, ngy
      do ix = 1, ngx
        np = grid(ix, iy)%np
        if (np == 0) cycle

        ! Extract particle data for this cell
        allocate(xs_cell(np), ys_cell(np), qs_cell(np))
        xs_cell = xs(grid(ix, iy)%plist)
        ys_cell = ys(grid(ix, iy)%plist)
        qs_cell = qs(grid(ix, iy)%plist)

        ! Compute multipole moments
        call compute_multipole_2d(np, xs_cell, ys_cell, qs_cell, &
                                 grid(ix, iy)%xc, grid(ix, iy)%yc, &
                                 p_order, grid(ix, iy)%M)

        deallocate(xs_cell, ys_cell, qs_cell)
      end do
    end do

  end subroutine p2m_step

  !----------------------------------------------------------------------------
  ! M2L step: Translate multipole to local for far-field interactions
  !----------------------------------------------------------------------------
  subroutine m2l_step(grid, ngx, ngy, p_order, near_range)
    type(grid_cell_2d), intent(inout) :: grid(:,:)
    integer, intent(in) :: ngx, ngy, p_order, near_range

    integer :: ix, iy, jx, jy
    real(dp) :: dx, dy, sep
    logical :: is_near

    ! For each target cell
    do iy = 1, ngy
      do ix = 1, ngx

        ! For each source cell
        do jy = 1, ngy
          do jx = 1, ngx
            if (ix == jx .and. iy == jy) cycle

            ! Check if in near field
            is_near = (abs(ix - jx) <= near_range) .and. &
                     (abs(iy - jy) <= near_range)

            if (is_near) cycle  ! Handle in P2P step

            ! Skip empty source cells
            if (grid(jx, jy)%np == 0) cycle

            ! M2L translation
            call m2l_translate_2d(grid(jx, jy)%M, &
                                 grid(jx, jy)%xc, grid(jx, jy)%yc, &
                                 grid(ix, iy)%xc, grid(ix, iy)%yc, &
                                 p_order, grid(ix, iy)%L)
          end do
        end do

      end do
    end do

  end subroutine m2l_step

  !----------------------------------------------------------------------------
  ! Evaluation step: L2P (far field) + P2P (near field)
  !----------------------------------------------------------------------------
  subroutine evaluate_step(grid, ngx, ngy, N, xs, ys, qs, potential, &
                          p_order, near_range, phi, fx, fy)
    type(grid_cell_2d), intent(in) :: grid(:,:)
    integer, intent(in) :: ngx, ngy, N, p_order, near_range
    real(dp), intent(in) :: xs(N), ys(N), qs(N)
    class(potential_type), intent(in) :: potential
    real(dp), intent(inout) :: phi(N), fx(N), fy(N)

    integer :: ix, iy, i, ip
    real(dp) :: phi_local, fx_local, fy_local

    ! Process each cell
    do iy = 1, ngy
      do ix = 1, ngx
        if (grid(ix, iy)%np == 0) cycle

        ! L2P: Evaluate local expansion at each particle
        do i = 1, grid(ix, iy)%np
          ip = grid(ix, iy)%plist(i)

          call evaluate_local_2d(grid(ix, iy)%L, &
                                grid(ix, iy)%xc, grid(ix, iy)%yc, &
                                xs(ip), ys(ip), p_order, &
                                phi_local, fx_local, fy_local)

          phi(ip) = phi(ip) + phi_local
          fx(ip) = fx(ip) + fx_local
          fy(ip) = fy(ip) + fy_local
        end do

        ! P2P: Direct calculation with near-field cells
        call p2p_near_field(grid, ix, iy, ngx, ngy, near_range, &
                           N, xs, ys, qs, potential, phi, fx, fy)
      end do
    end do

  end subroutine evaluate_step

  !----------------------------------------------------------------------------
  ! P2P: Direct near-field interactions
  !----------------------------------------------------------------------------
  subroutine p2p_near_field(grid, ix, iy, ngx, ngy, near_range, &
                           N, xs, ys, qs, potential, phi, fx, fy)
    type(grid_cell_2d), intent(in) :: grid(:,:)
    integer, intent(in) :: ix, iy, ngx, ngy, near_range, N
    real(dp), intent(in) :: xs(N), ys(N), qs(N)
    class(potential_type), intent(in) :: potential
    real(dp), intent(inout) :: phi(N), fx(N), fy(N)

    integer :: jx, jy, i, j, ip, jp
    real(dp) :: rx, ry, r, fx_ij, fy_ij, phi_ij
    logical :: same_cell

    ! Loop over nearby cells
    do jy = max(1, iy - near_range), min(ngy, iy + near_range)
      do jx = max(1, ix - near_range), min(ngx, ix + near_range)
        if (grid(jx, jy)%np == 0) cycle

        same_cell = (ix == jx) .and. (iy == jy)

        ! Particle-particle interactions
        do i = 1, grid(ix, iy)%np
          ip = grid(ix, iy)%plist(i)

          do j = 1, grid(jx, jy)%np
            jp = grid(jx, jy)%plist(j)

            ! Skip self-interaction
            if (same_cell .and. ip == jp) cycle

            ! Compute separation
            rx = xs(ip) - xs(jp)
            ry = ys(ip) - ys(jp)
            r = sqrt(rx*rx + ry*ry)

            if (r < TINY) cycle

            ! Compute potential and force using potential interface
            phi_ij = potential%value(r, qs(jp))
            call potential%force(rx, ry, qs(jp), fx_ij, fy_ij)

            phi(ip) = phi(ip) + phi_ij
            fx(ip) = fx(ip) + fx_ij
            fy(ip) = fy(ip) + fy_ij
          end do
        end do

      end do
    end do

  end subroutine p2p_near_field

  !----------------------------------------------------------------------------
  ! Cleanup grid
  !----------------------------------------------------------------------------
  subroutine cleanup_grid(grid, ngx, ngy)
    type(grid_cell_2d), intent(inout) :: grid(:,:)
    integer, intent(in) :: ngx, ngy

    integer :: ix, iy

    do iy = 1, ngy
      do ix = 1, ngx
        if (allocated(grid(ix, iy)%plist)) deallocate(grid(ix, iy)%plist)
        if (allocated(grid(ix, iy)%M)) deallocate(grid(ix, iy)%M)
        if (allocated(grid(ix, iy)%L)) deallocate(grid(ix, iy)%L)
      end do
    end do

  end subroutine cleanup_grid

end module grid_fmm_mod
