module cartesian_multipole_mod
  !==============================================================================
  ! Cartesian multipole expansion for 2D problems
  !
  ! For 2D Coulomb potential 1/r, we use Cartesian multipole expansion:
  !
  ! Multipole expansion (source -> multipole):
  !   M_{n,m} = sum_i q_i * r_i^n * T_n^m(theta_i)
  !
  ! Local expansion (multipole -> field point):
  !   phi = sum_{n,m} L_{n,m} / r^{n+1} * T_n^m(theta)
  !
  ! Where T_n^m are Chebyshev-like polynomials:
  !   T_n^m(theta) = cos(m*theta) for m >= 0
  !   T_n^m(theta) = sin(|m|*theta) for m < 0
  !
  ! In Cartesian coordinates:
  !   x = r * cos(theta)
  !   y = r * sin(theta)
  !
  ! Using complex notation: z = x + iy = r * e^(i*theta)
  ! Then: z^n = r^n * e^(i*n*theta) = r^n * (cos(n*theta) + i*sin(n*theta))
  !==============================================================================
  use types_mod
  implicit none
  private

  public :: compute_multipole_2d
  public :: compute_local_2d
  public :: m2m_translate_2d
  public :: m2l_translate_2d
  public :: l2l_translate_2d
  public :: evaluate_local_2d
  public :: precompute_factorials

  ! Precomputed factorials for normalization
  real(dp), allocatable, save :: factorials(:)
  integer, save :: max_factorial = -1

contains

  !------------------------------------------------------------------------------
  ! Precompute factorials up to n
  !------------------------------------------------------------------------------
  subroutine precompute_factorials(n)
    integer, intent(in) :: n
    integer :: i

    if (allocated(factorials)) then
      if (max_factorial >= n) return
      deallocate(factorials)
    end if

    max_factorial = n
    allocate(factorials(0:n))

    factorials(0) = 1.0_dp
    do i = 1, n
      factorials(i) = factorials(i-1) * real(i, dp)
    end do
  end subroutine precompute_factorials

  !------------------------------------------------------------------------------
  ! Compute binomial coefficient C(n, k)
  !------------------------------------------------------------------------------
  pure recursive function binomial(n, k) result(c)
    integer, intent(in) :: n, k
    real(dp) :: c
    integer :: i

    if (k < 0 .or. k > n) then
      c = 0.0_dp
      return
    end if

    if (k == 0 .or. k == n) then
      c = 1.0_dp
      return
    end if

    ! Use smaller range
    if (k > n - k) then
      c = binomial(n, n - k)
      return
    end if

    c = 1.0_dp
    do i = 0, k - 1
      c = c * real(n - i, dp) / real(i + 1, dp)
    end do
  end function binomial

  !------------------------------------------------------------------------------
  ! Compute P2M: Particle to Multipole
  ! M_{n,m} = sum_i q_i * (x_i - xc)^n * (y_i - yc)^m
  !------------------------------------------------------------------------------
  subroutine compute_multipole_2d(np, xs, ys, qs, xc, yc, p_order, Mout)
    integer, intent(in) :: np, p_order
    real(dp), intent(in) :: xs(np), ys(np), qs(np)
    real(dp), intent(in) :: xc, yc
    real(dp), intent(out) :: Mout(0:p_order, 0:p_order)

    integer :: i, n, m
    real(dp) :: dx, dy, dxn, dym

    ! Initialize moments to zero
    Mout = 0.0_dp

    ! Compute moments
    do i = 1, np
      dx = xs(i) - xc
      dy = ys(i) - yc

      ! Iterate over orders
      do n = 0, p_order
        dxn = dx**n
        do m = 0, p_order - n
          dym = dy**m
          Mout(n, m) = Mout(n, m) + qs(i) * dxn * dym
        end do
      end do
    end do
  end subroutine compute_multipole_2d

  !------------------------------------------------------------------------------
  ! M2M: Multipole to Multipole translation
  ! Translate multipole from child (xc, yc) to parent (xp, yp)
  !
  ! M_parent(n,m) = sum_{k=0}^n sum_{l=0}^m C(n,k)*C(m,l) * (dx)^(n-k) * (dy)^(m-l) * M_child(k,l)
  !------------------------------------------------------------------------------
  subroutine m2m_translate_2d(M_child, xc, yc, xp, yp, p_order, M_parent)
    real(dp), intent(in) :: M_child(0:p_order, 0:p_order)
    real(dp), intent(in) :: xc, yc, xp, yp
    integer, intent(in) :: p_order
    real(dp), intent(inout) :: M_parent(0:p_order, 0:p_order)

    integer :: n, m, k, l
    real(dp) :: dx, dy, dx_pow, dy_pow, coeff

    dx = xc - xp
    dy = yc - yp

    ! Translate moments using Taylor expansion
    do n = 0, p_order
      do m = 0, p_order - n
        do k = 0, n
          dx_pow = dx**(n - k)
          do l = 0, m
            dy_pow = dy**(m - l)
            coeff = binomial(n, k) * binomial(m, l) * dx_pow * dy_pow

            M_parent(n, m) = M_parent(n, m) + coeff * M_child(k, l)
          end do
        end do
      end do
    end do
  end subroutine m2m_translate_2d

  !------------------------------------------------------------------------------
  ! M2L: Multipole to Local translation
  ! Translate multipole at (xs, ys) to local expansion at (xt, yt)
  !
  ! For well-separated cells, convert multipole to local expansion
  !------------------------------------------------------------------------------
  subroutine m2l_translate_2d(M_source, xs, ys, xt, yt, p_order, L_target)
    real(dp), intent(in) :: M_source(0:p_order, 0:p_order)
    real(dp), intent(in) :: xs, ys, xt, yt
    integer, intent(in) :: p_order
    real(dp), intent(inout) :: L_target(0:p_order, 0:p_order)

    integer :: n, m, k, l, j
    real(dp) :: dx, dy, r, r2, r_inv
    real(dp) :: dx_n, dy_m, r_pow
    complex(dp) :: z, z_pow, z_k, contrib
    complex(dp) :: M_complex, L_complex

    dx = xs - xt
    dy = ys - yt
    r2 = dx*dx + dy*dy
    r = sqrt(r2)

    if (r < SMALL) return

    r_inv = 1.0_dp / r

    ! Use complex arithmetic for efficiency
    ! z = dx + i*dy, z^n = r^n * e^(i*n*theta)
    z = cmplx(dx, dy, kind=dp)

    ! Compute L_{j} from M_{k} using:
    ! L_j = sum_k M_k * (-1)^k / (z^(j+k+1))
    do j = 0, p_order
      do n = 0, j
        do m = 0, j - n
          L_complex = cmplx(0.0_dp, 0.0_dp, kind=dp)

          do k = 0, p_order
            z_k = z**(j + k + 1)
            z_pow = 1.0_dp / z_k

            do l = 0, p_order - k
              ! Convert (k, l) Cartesian to complex
              M_complex = cmplx(M_source(k, l), 0.0_dp, kind=dp) * &
                          cmplx(real(k + l, dp), 0.0_dp, kind=dp)

              ! Compute contribution using convolution
              contrib = M_complex * z_pow * (-1.0_dp)**k / factorials(k + l)

              L_complex = L_complex + contrib
            end do
          end do

          ! Extract real part
          L_target(n, m) = L_target(n, m) + real(L_complex, kind=dp)
        end do
      end do
    end do
  end subroutine m2l_translate_2d

  !------------------------------------------------------------------------------
  ! L2L: Local to Local translation
  ! Translate local expansion from parent (xp, yp) to child (xc, yc)
  !
  ! L_child(n,m) = sum_{k=n}^p sum_{l=m}^{p-k} C(k,n)*C(l,m) * (dx)^(k-n) * (dy)^(l-m) * L_parent(k,l)
  !------------------------------------------------------------------------------
  subroutine l2l_translate_2d(L_parent, xp, yp, xc, yc, p_order, L_child)
    real(dp), intent(in) :: L_parent(0:p_order, 0:p_order)
    real(dp), intent(in) :: xp, yp, xc, yc
    integer, intent(in) :: p_order
    real(dp), intent(inout) :: L_child(0:p_order, 0:p_order)

    integer :: n, m, k, l
    real(dp) :: dx, dy, dx_pow, dy_pow, coeff

    dx = xc - xp
    dy = yc - yp

    ! Translate local expansion using Taylor shift
    do n = 0, p_order
      do m = 0, p_order - n
        do k = n, p_order
          dx_pow = dx**(k - n)
          do l = m, p_order - k
            dy_pow = dy**(l - m)
            coeff = binomial(k, n) * binomial(l, m) * dx_pow * dy_pow

            L_child(n, m) = L_child(n, m) + coeff * L_parent(k, l)
          end do
        end do
      end do
    end do
  end subroutine l2l_translate_2d

  !------------------------------------------------------------------------------
  ! Initialize local expansion from multipole moments (for far field)
  ! This converts multipole M at (xm, ym) to local L at (xl, yl)
  !------------------------------------------------------------------------------
  subroutine compute_local_2d(M, xm, ym, xl, yl, p_order, L)
    real(dp), intent(in) :: M(0:p_order, 0:p_order)
    real(dp), intent(in) :: xm, ym, xl, yl
    integer, intent(in) :: p_order
    real(dp), intent(out) :: L(0:p_order, 0:p_order)

    L = 0.0_dp
    call m2l_translate_2d(M, xm, ym, xl, yl, p_order, L)
  end subroutine compute_local_2d

  !------------------------------------------------------------------------------
  ! Evaluate local expansion at point (x, y) relative to center (xc, yc)
  ! Returns potential phi and force (fx, fy)
  !------------------------------------------------------------------------------
  subroutine evaluate_local_2d(L, xc, yc, x, y, p_order, phi, fx, fy)
    real(dp), intent(in) :: L(0:p_order, 0:p_order)
    real(dp), intent(in) :: xc, yc, x, y
    integer, intent(in) :: p_order
    real(dp), intent(out) :: phi, fx, fy

    integer :: n, m
    real(dp) :: dx, dy, dx_n, dy_m, dx_n1, dy_m1

    dx = x - xc
    dy = y - yc

    phi = 0.0_dp
    fx = 0.0_dp
    fy = 0.0_dp

    ! Evaluate phi = sum L(n,m) * dx^n * dy^m
    do n = 0, p_order
      dx_n = dx**n
      if (n > 0) then
        dx_n1 = real(n, dp) * dx**(n-1)
      else
        dx_n1 = 0.0_dp
      end if

      do m = 0, p_order - n
        dy_m = dy**m
        if (m > 0) then
          dy_m1 = real(m, dp) * dy**(m-1)
        else
          dy_m1 = 0.0_dp
        end if

        phi = phi + L(n, m) * dx_n * dy_m

        ! F = -grad(phi)
        fx = fx - L(n, m) * dx_n1 * dy_m
        fy = fy - L(n, m) * dx_n * dy_m1
      end do
    end do
  end subroutine evaluate_local_2d

end module cartesian_multipole_mod
