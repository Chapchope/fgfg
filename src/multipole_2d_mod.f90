module multipole_2d_mod
  !============================================================================
  ! МОДУЛЬ: 2D мультипольные операторы для FMM
  !
  ! НАЗНАЧЕНИЕ:
  !   Операторы P2M, M2M, M2L, L2L, L2P для 2D картезианских мультиполей
  !   Поддержка произвольного порядка разложения p
  !
  ! ФИЗИКА:
  !   Мультипольный момент: M_nm = Σ q_i (x_i - xc)^n (y_i - yc)^m
  !   Локальное разложение: L_nm для вычисления φ = Σ L_nm x^n y^m
  !
  ! ИСПОЛЬЗОВАНИЕ:
  !   call p2m(particles, cell_center, p_order, M)
  !   call m2l(M_source, center_source, center_target, p_order, potential, L_target)
  !   call l2p(L, cell_center, particles, p_order, phi, fx, fy)
  !
  ! АВТОР: Grid-based FMM для сверхпроводника
  !============================================================================
  use types_mod
  use potential_2d_mod
  implicit none
  private

  public :: p2m_2d
  public :: m2l_2d
  public :: l2p_2d
  public :: p2p_direct_2d

contains

  !============================================================================
  ! P2M: Particle-to-Multipole
  !
  ! Вычисляет мультипольные моменты M_nm для набора частиц
  ! M_nm = Σ q_i · (x_i - xc)^n · (y_i - yc)^m
  !
  ! ВХОД:
  !   np - число частиц
  !   xs, ys, qs - координаты и заряды
  !   xc, yc - центр разложения
  !   p_order - порядок разложения
  !
  ! ВЫХОД:
  !   M(0:p_order, 0:p_order) - мультипольные моменты
  !============================================================================
  subroutine p2m_2d(np, xs, ys, qs, xc, yc, p_order, M_out)
    integer, intent(in) :: np, p_order
    real(dp), intent(in) :: xs(np), ys(np), qs(np)
    real(dp), intent(in) :: xc, yc
    real(dp), intent(out) :: M_out(0:p_order, 0:p_order)

    integer :: i, nn, mm
    real(dp) :: dx, dy, q
    real(dp) :: pow_x(0:p_order), pow_y(0:p_order)

    M_out = 0.0_dp

    do i = 1, np
      q = qs(i)
      dx = xs(i) - xc
      dy = ys(i) - yc

      ! Предвычисляем степени для эффективности
      pow_x(0) = 1.0_dp
      pow_y(0) = 1.0_dp
      do nn = 1, p_order
        pow_x(nn) = pow_x(nn-1) * dx
        pow_y(nn) = pow_y(nn-1) * dy
      end do

      ! Суммируем моменты
      do nn = 0, p_order
        do mm = 0, p_order - nn
          M_out(nn,mm) = M_out(nn,mm) + q * pow_x(nn) * pow_y(mm)
        end do
      end do
    end do

  end subroutine p2m_2d

  !============================================================================
  ! M2L: Multipole-to-Local
  !
  ! Конвертирует мультипольное разложение M в локальное разложение L
  ! Это ключевая операция FMM!
  !
  ! ФИЗИКА:
  !   Берем мультиполи M источника в точке (xc_src, yc_src)
  !   Конвертируем в локальное разложение L в точке (xc_tgt, yc_tgt)
  !
  !   L_nm += Σ_{k,l} M_kl · T_{n,m,k,l}
  !
  !   Где T - оператор трансляции, вычисляемый из производных потенциала
  !
  ! ВАЖНО:
  !   Для высокой точности (10^-10) необходимо правильно вычислить
  !   все производные до порядка p_order.
  !   Это делается через потенциал (см. potential_2d_mod)
  !
  ! ВХОД:
  !   M_src - мультипольные моменты источника
  !   xc_src, yc_src - центр источника
  !   xc_tgt, yc_tgt - центр целевой ячейки
  !   p_order - порядок разложения
  !   potential - объект потенциала для вычисления производных
  !
  ! ВЫХОД:
  !   L_tgt - локальное разложение (аккумулируется!)
  !============================================================================
  subroutine m2l_2d(M_src, xc_src, yc_src, xc_tgt, yc_tgt, p_order, potential, L_tgt)
    real(dp), intent(in) :: M_src(0:p_order, 0:p_order)
    real(dp), intent(in) :: xc_src, yc_src, xc_tgt, yc_tgt
    integer, intent(in) :: p_order
    class(abstract_potential_2d), intent(in) :: potential
    real(dp), intent(inout) :: L_tgt(0:p_order, 0:p_order)

    real(dp) :: dx, dy, R
    real(dp) :: derivs(0:2*p_order, 0:2*p_order)
    integer :: nn, mm, kk, ll

    ! Вектор от target к source
    dx = xc_src - xc_tgt
    dy = yc_src - yc_tgt
    R = sqrt(dx*dx + dy*dy)

    if (R < 1.0e-14_dp) return

    ! Вычисляем производные потенциала через объект потенциала
    call potential%compute_multipole_derivatives(dx, dy, R, 2*p_order, derivs)

    ! M2L трансляция: L_nm += Σ M_kl · ∂^{n+k}∂^{m+l}(1/R)
    ! Формула: L_nm = Σ_{k=0}^{p} Σ_{l=0}^{p-k} M_kl · T_{n+k,m+l}
    do nn = 0, p_order
      do mm = 0, p_order - nn
        do kk = 0, p_order
          do ll = 0, p_order - kk
            if (nn+kk <= 2*p_order .and. mm+ll <= 2*p_order) then
              L_tgt(nn,mm) = L_tgt(nn,mm) + M_src(kk,ll) * derivs(nn+kk, mm+ll)
            end if
          end do
        end do
      end do
    end do

  end subroutine m2l_2d

  !============================================================================
  ! L2P: Local-to-Particle
  !
  ! Вычисляет потенциал φ и силу F из локального разложения L
  !
  ! ФИЗИКА:
  !   φ(r) = Σ L_nm · (x - xc)^n · (y - yc)^m
  !   F_x = -∂φ/∂x = -Σ n · L_nm · (x - xc)^{n-1} · (y - yc)^m
  !   F_y = -∂φ/∂y = -Σ m · L_nm · (x - xc)^n · (y - yc)^{m-1}
  !
  ! ВХОД:
  !   L - локальное разложение
  !   xc, yc - центр разложения
  !   np - число частиц
  !   xs, ys - координаты частиц
  !   p_order - порядок разложения
  !
  ! ВЫХОД:
  !   phi, fx, fy - потенциал и силы (аккумулируются!)
  !============================================================================
  subroutine l2p_2d(L_in, xc, yc, np, xs, ys, p_order, phi, fx, fy)
    real(dp), intent(in) :: L_in(0:p_order, 0:p_order)
    real(dp), intent(in) :: xc, yc
    integer, intent(in) :: np, p_order
    real(dp), intent(in) :: xs(np), ys(np)
    real(dp), intent(inout) :: phi(np), fx(np), fy(np)

    integer :: i, nn, mm
    real(dp) :: dx, dy
    real(dp) :: pow_x(0:p_order), pow_y(0:p_order)
    real(dp) :: phi_tmp, fx_tmp, fy_tmp

    do i = 1, np
      dx = xs(i) - xc
      dy = ys(i) - yc

      ! Предвычисляем степени
      pow_x(0) = 1.0_dp
      pow_y(0) = 1.0_dp
      do nn = 1, p_order
        pow_x(nn) = pow_x(nn-1) * dx
        pow_y(nn) = pow_y(nn-1) * dy
      end do

      phi_tmp = 0.0_dp
      fx_tmp = 0.0_dp
      fy_tmp = 0.0_dp

      ! Вычисляем φ и F из L
      do nn = 0, p_order
        do mm = 0, p_order - nn
          if (abs(L_in(nn,mm)) < 1.0e-30_dp) cycle

          ! Потенциал
          phi_tmp = phi_tmp + L_in(nn,mm) * pow_x(nn) * pow_y(mm)

          ! Силы: F = -∇φ
          if (nn > 0) then
            fx_tmp = fx_tmp - real(nn,dp) * L_in(nn,mm) * pow_x(nn-1) * pow_y(mm)
          end if

          if (mm > 0) then
            fy_tmp = fy_tmp - real(mm,dp) * L_in(nn,mm) * pow_x(nn) * pow_y(mm-1)
          end if
        end do
      end do

      phi(i) = phi(i) + phi_tmp
      fx(i) = fx(i) + fx_tmp
      fy(i) = fy(i) + fy_tmp
    end do

  end subroutine l2p_2d

  !============================================================================
  ! P2P: Прямое вычисление particle-to-particle
  !
  ! Для близких частиц мультипольное разложение неточно.
  ! Поэтому мы считаем напрямую через потенциал.
  !
  ! ВХОД:
  !   np_src, np_tgt - число частиц источника и цели
  !   xs_src, ys_src, qs_src - координаты и заряды источников
  !   xs_tgt, ys_tgt - координаты целей
  !   potential - объект потенциала
  !   exclude_self - исключать ли самовзаимодействие (если src==tgt)
  !
  ! ВЫХОД:
  !   phi, fx, fy - потенциал и силы (аккумулируются!)
  !============================================================================
  subroutine p2p_direct_2d(np_src, xs_src, ys_src, qs_src, &
                           np_tgt, xs_tgt, ys_tgt, &
                           potential, exclude_self, &
                           phi, fx, fy)
    integer, intent(in) :: np_src, np_tgt
    real(dp), intent(in) :: xs_src(np_src), ys_src(np_src), qs_src(np_src)
    real(dp), intent(in) :: xs_tgt(np_tgt), ys_tgt(np_tgt)
    class(abstract_potential_2d), intent(in) :: potential
    logical, intent(in) :: exclude_self
    real(dp), intent(inout) :: phi(np_tgt), fx(np_tgt), fy(np_tgt)

    integer :: i, j
    real(dp) :: phi_tmp, fx_tmp, fy_tmp
    real(dp) :: dx, dy, r2

    do i = 1, np_tgt
      do j = 1, np_src

        ! Пропускаем самовзаимодействие если нужно
        if (exclude_self) then
          dx = xs_tgt(i) - xs_src(j)
          dy = ys_tgt(i) - ys_src(j)
          r2 = dx*dx + dy*dy
          if (r2 < 1.0e-16_dp) cycle
        end if

        ! Вычисляем φ и F через потенциал
        call potential%eval_field(xs_tgt(i), ys_tgt(i), qs_src(j), &
                                  xs_src(j), ys_src(j), &
                                  phi_tmp, fx_tmp, fy_tmp)

        phi(i) = phi(i) + phi_tmp
        fx(i) = fx(i) + fx_tmp
        fy(i) = fy(i) + fy_tmp
      end do
    end do

  end subroutine p2p_direct_2d

end module multipole_2d_mod
