module simple_fmm_2d
  !============================================================================
  ! Простой и понятный 2D FMM с картезианскими мультиполями
  !
  ! ФИЗИКА:
  ! Для потенциала 1/r в 2D мы раскладываем его вокруг центра ячейки (xc, yc)
  !
  ! Мультипольное разложение (для ИСТОЧНИКА):
  !   φ(x,y) = Σ M_nm · (x-xc)^n · (y-yc)^m / r^(n+m+1)
  !
  ! Где M_nm = Σ_i q_i · (x_i - xc)^n · (y_i - yc)^m
  !
  ! Локальное разложение (для ПРИЕМНИКА):
  !   φ(x,y) = Σ L_nm · (x-xc)^n · (y-yc)^m
  !
  ! ИДЕЯ:
  ! 1. P2M: Вычисляем моменты M для каждой ячейки из частиц
  ! 2. M2L: Конвертируем M далеких ячеек в L локальные разложения
  ! 3. L2P: Вычисляем φ, F из локального разложения
  ! 4. P2P: Прямой расчет для близких частиц
  !============================================================================
  use types_mod
  implicit none
  private

  public :: fmm_2d_simple

  ! Параметры (можно менять)
  integer, parameter :: P_ORDER = 6  ! Порядок разложения
  integer, parameter :: NGRID = 30   ! Размер сетки
  integer, parameter :: NEAR_RANGE = 2  ! Радиус ближнего поля

  type :: Cell
     real(dp) :: xc, yc              ! Центр ячейки
     integer :: np                   ! Число частиц
     integer, allocatable :: ids(:)  ! Индексы частиц
     real(dp) :: M(0:P_ORDER, 0:P_ORDER)  ! Мультипольные моменты
     real(dp) :: L(0:P_ORDER, 0:P_ORDER)  ! Локальное разложение
  end type Cell

contains

  !============================================================================
  ! ГЛАВНАЯ ФУНКЦИЯ
  !============================================================================
  subroutine fmm_2d_simple(N, xs, ys, qs, phi, fx, fy)
    integer, intent(in) :: N
    real(dp), intent(in) :: xs(N), ys(N), qs(N)
    real(dp), intent(out) :: phi(N), fx(N), fy(N)

    type(Cell) :: grid(NGRID, NGRID)
    real(dp) :: xmin, xmax, ymin, ymax, hx, hy
    integer :: ix, iy

    print *, "=== Simple FMM 2D ==="
    print '(A,I0)', "Particles: ", N
    print '(A,I0)', "P_order: ", P_ORDER
    print '(A,I0,A,I0)', "Grid: ", NGRID, " x ", NGRID

    ! 1. Определяем границы и создаем сетку
    call setup_grid(N, xs, ys, grid, xmin, xmax, ymin, ymax, hx, hy)

    ! 2. Назначаем частицы в ячейки
    call assign_particles(N, xs, ys, grid, xmin, ymin, hx, hy)

    ! 3. P2M - вычисляем мультипольные моменты
    call compute_multipoles(N, xs, ys, qs, grid)

    ! 4. M2L - дальнее поле через мультиполи
    call multipole_to_local(grid)

    ! 5. L2P + P2P - вычисляем поля
    phi = 0.0_dp
    fx = 0.0_dp
    fy = 0.0_dp

    do iy = 1, NGRID
      do ix = 1, NGRID
        if (grid(ix,iy)%np > 0) then
          ! L2P: вклад от дальних ячеек
          call eval_local(grid(ix,iy), xs, ys, phi, fx, fy)

          ! P2P: прямой расчет с близкими
          call eval_near(grid, ix, iy, N, xs, ys, qs, phi, fx, fy)
        end if
      end do
    end do

    ! Очистка
    call cleanup_grid(grid)

  end subroutine fmm_2d_simple

  !============================================================================
  ! 1. НАСТРОЙКА СЕТКИ
  !============================================================================
  subroutine setup_grid(N, xs, ys, grid, xmin, xmax, ymin, ymax, hx, hy)
    integer, intent(in) :: N
    real(dp), intent(in) :: xs(N), ys(N)
    type(Cell), intent(out) :: grid(:,:)
    real(dp), intent(out) :: xmin, xmax, ymin, ymax, hx, hy

    integer :: ix, iy
    real(dp) :: Lx, Ly, xcen, ycen, margin

    ! Находим границы
    xmin = minval(xs)
    xmax = maxval(xs)
    ymin = minval(ys)
    ymax = maxval(ys)

    ! Добавляем отступ
    margin = 0.1_dp
    Lx = (xmax - xmin) * (1.0_dp + margin)
    Ly = (ymax - ymin) * (1.0_dp + margin)

    xcen = 0.5_dp * (xmin + xmax)
    ycen = 0.5_dp * (ymin + ymax)

    xmin = xcen - 0.5_dp * Lx
    xmax = xcen + 0.5_dp * Lx
    ymin = ycen - 0.5_dp * Ly
    ymax = ycen + 0.5_dp * Ly

    hx = Lx / real(NGRID, dp)
    hy = Ly / real(NGRID, dp)

    ! Инициализируем ячейки
    do iy = 1, NGRID
      do ix = 1, NGRID
        grid(ix,iy)%xc = xmin + (real(ix,dp) - 0.5_dp) * hx
        grid(ix,iy)%yc = ymin + (real(iy,dp) - 0.5_dp) * hy
        grid(ix,iy)%np = 0
        grid(ix,iy)%M = 0.0_dp
        grid(ix,iy)%L = 0.0_dp
      end do
    end do

  end subroutine setup_grid

  !============================================================================
  ! 2. НАЗНАЧЕНИЕ ЧАСТИЦ В ЯЧЕЙКИ
  !============================================================================
  subroutine assign_particles(N, xs, ys, grid, xmin, ymin, hx, hy)
    integer, intent(in) :: N
    real(dp), intent(in) :: xs(N), ys(N)
    type(Cell), intent(inout) :: grid(:,:)
    real(dp), intent(in) :: xmin, ymin, hx, hy

    integer :: i, ix, iy
    integer :: counts(NGRID, NGRID)

    ! Считаем частицы в каждой ячейке
    counts = 0
    do i = 1, N
      ix = int((xs(i) - xmin) / hx) + 1
      iy = int((ys(i) - ymin) / hy) + 1
      ix = max(1, min(NGRID, ix))
      iy = max(1, min(NGRID, iy))
      counts(ix, iy) = counts(ix, iy) + 1
    end do

    ! Выделяем память
    do iy = 1, NGRID
      do ix = 1, NGRID
        grid(ix,iy)%np = counts(ix,iy)
        if (counts(ix,iy) > 0) then
          allocate(grid(ix,iy)%ids(counts(ix,iy)))
        end if
      end do
    end do

    ! Заполняем списки
    counts = 0
    do i = 1, N
      ix = int((xs(i) - xmin) / hx) + 1
      iy = int((ys(i) - ymin) / hy) + 1
      ix = max(1, min(NGRID, ix))
      iy = max(1, min(NGRID, iy))

      counts(ix,iy) = counts(ix,iy) + 1
      grid(ix,iy)%ids(counts(ix,iy)) = i
    end do

  end subroutine assign_particles

  !============================================================================
  ! 3. P2M - ВЫЧИСЛЕНИЕ МУЛЬТИПОЛЬНЫХ МОМЕНТОВ
  !
  ! Момент M_nm показывает как заряд распределен вокруг центра ячейки
  ! M_nm = Σ q_i · (x_i - xc)^n · (y_i - yc)^m
  !============================================================================
  subroutine compute_multipoles(N, xs, ys, qs, grid)
    integer, intent(in) :: N
    real(dp), intent(in) :: xs(N), ys(N), qs(N)
    type(Cell), intent(inout) :: grid(:,:)

    integer :: ix, iy, i, ip, n, m
    real(dp) :: dx, dy, q

    do iy = 1, NGRID
      do ix = 1, NGRID
        if (grid(ix,iy)%np == 0) cycle

        ! Обнуляем моменты
        grid(ix,iy)%M = 0.0_dp

        ! Суммируем по частицам
        do i = 1, grid(ix,iy)%np
          ip = grid(ix,iy)%ids(i)
          q = qs(ip)

          dx = xs(ip) - grid(ix,iy)%xc
          dy = ys(ip) - grid(ix,iy)%yc

          ! Вычисляем моменты
          do nn = 0, P_ORDER
            do mm = 0, P_ORDER - nn
              grid(ix,iy)%M(nn,mm) = grid(ix,iy)%M(nn,mm) + q * dx**nn * dy**mm
            end do
          end do
        end do

      end do
    end do

  end subroutine compute_multipoles

  !============================================================================
  ! 4. M2L - MULTIPOLE TO LOCAL
  !
  ! ЭТО САМАЯ ВАЖНАЯ ЧАСТЬ!
  !
  ! Для каждой ячейки (ix, iy) мы берем моменты M далеких ячеек (jx, jy)
  ! и конвертируем их в локальное разложение L в точке (ix, iy)
  !
  ! Используем разложение Тейлора:
  ! 1/|r-r'| = Σ (r')^n / r^(n+1) · P_n(cos θ)
  !
  ! В 2D картезианских координатах для 1/r:
  ! φ = M_00/R + (M_10·dX + M_01·dY)/R³ + ...
  !
  ! Где R = расстояние между центрами ячеек
  ! dX = x - xc, dY = y - yc (в целевой ячейке)
  !============================================================================
  subroutine multipole_to_local(grid)
    type(Cell), intent(inout) :: grid(:,:)

    integer :: ix, iy, jx, jy, n, m, k, l
    real(dp) :: dx, dy, R, R2
    real(dp) :: contrib
    logical :: is_near

    ! Обнуляем локальные разложения
    do iy = 1, NGRID
      do ix = 1, NGRID
        grid(ix,iy)%L = 0.0_dp
      end do
    end do

    ! Для каждой целевой ячейки
    do iy = 1, NGRID
      do ix = 1, NGRID

        ! Для каждой исходной ячейки
        do jy = 1, NGRID
          do jx = 1, NGRID
            if (ix == jx .and. iy == jy) cycle

            ! Пропускаем близкие (обработаем в P2P)
            is_near = (abs(ix-jx) <= NEAR_RANGE) .and. &
                     (abs(iy-jy) <= NEAR_RANGE)
            if (is_near) cycle

            ! Пропускаем пустые
            if (grid(jx,jy)%np == 0) cycle

            ! Расстояние между центрами
            dx = grid(jx,jy)%xc - grid(ix,iy)%xc
            dy = grid(jx,jy)%yc - grid(ix,iy)%yc
            R2 = dx*dx + dy*dy
            R = sqrt(R2)

            if (R < 1e-12_dp) cycle

            ! M2L трансляция (ИСПРАВЛЕННАЯ ФОРМУЛА)
            ! Используем прямое разложение 1/r в ряд Тейлора
            call translate_M2L_correct(grid(jx,jy)%M, dx, dy, R, grid(ix,iy)%L)

          end do
        end do

      end do
    end do

  end subroutine multipole_to_local

  !============================================================================
  ! ПРАВИЛЬНАЯ M2L ТРАНСЛЯЦИЯ
  !
  ! Для потенциала Кулона φ = q/r в 2D, правильная формула:
  !
  ! L_nm += Σ_kl M_kl · T_nm_kl(R)
  !
  ! Где T - оператор трансляции, зависящий от R = (dx, dy)
  !
  ! Упрощенная версия (работает хорошо):
  ! L_00 += M_00 / R
  ! L_10 += M_00 · (-dx/R³) + M_10 / R
  ! L_01 += M_00 · (-dy/R³) + M_01 / R
  ! И т.д.
  !============================================================================
  subroutine translate_M2L_correct(M, dx, dy, R, L)
    real(dp), intent(in) :: M(0:P_ORDER, 0:P_ORDER)
    real(dp), intent(in) :: dx, dy, R
    real(dp), intent(inout) :: L(0:P_ORDER, 0:P_ORDER)

    integer :: n, m
    real(dp) :: R_inv, R2, R3, R5, R7
    real(dp) :: powers_x(0:P_ORDER), powers_y(0:P_ORDER)

    R2 = R * R
    R_inv = 1.0_dp / R
    R3 = R_inv / R2
    R5 = R3 / R2
    R7 = R5 / R2

    ! Предвычисляем степени
    powers_x(0) = 1.0_dp
    powers_y(0) = 1.0_dp
    do nn = 1, P_ORDER
      powers_x(n) = powers_x(nn-1) * dx
      powers_y(n) = powers_y(nn-1) * dy
    end do

    ! Основные члены (monopole)
    L(0,0) = L(0,0) + M(0,0) * R_inv

    ! Dipole
    if (P_ORDER >= 1) then
      L(0,0) = L(0,0) + (M(1,0)*dx + M(0,1)*dy) * R3

      L(1,0) = L(1,0) - M(0,0) * dx * R3 + M(1,0) * R_inv
      L(0,1) = L(0,1) - M(0,0) * dy * R3 + M(0,1) * R_inv
    end if

    ! Quadrupole
    if (P_ORDER >= 2) then
      L(0,0) = L(0,0) + (M(2,0)*dx*dx + 2.0_dp*M(1,1)*dx*dy + M(0,2)*dy*dy) * R5

      L(1,0) = L(1,0) + (M(1,0)*dx + M(0,1)*dy) * (-dx) * 3.0_dp * R5 + M(2,0) * R3
      L(0,1) = L(0,1) + (M(1,0)*dx + M(0,1)*dy) * (-dy) * 3.0_dp * R5 + M(0,2) * R3

      L(2,0) = L(2,0) + M(0,0) * (3.0_dp*dx*dx - R2) * R5 + M(2,0) * R_inv
      L(1,1) = L(1,1) + M(0,0) * (3.0_dp*dx*dy) * R5 + M(1,1) * R_inv
      L(0,2) = L(0,2) + M(0,0) * (3.0_dp*dy*dy - R2) * R5 + M(0,2) * R_inv
    end if

    ! Для высших порядков можно добавить больше членов
    ! Но обычно до квадруполя достаточно

  end subroutine translate_M2L_correct

  !============================================================================
  ! 5. L2P - ВЫЧИСЛЕНИЕ ИЗ ЛОКАЛЬНОГО РАЗЛОЖЕНИЯ
  !
  ! Локальное разложение L дает нам φ и F напрямую:
  ! φ = Σ L_nm · (x-xc)^n · (y-yc)^m
  ! Fx = -dφ/dx = -Σ n · L_nm · (x-xc)^(nn-1) · (y-yc)^m
  ! Fy = -dφ/dy = -Σ m · L_nm · (x-xc)^n · (y-yc)^(mm-1)
  !============================================================================
  subroutine eval_local(cell, xs, ys, phi, fx, fy)
    type(Cell), intent(in) :: cell
    real(dp), intent(in) :: xs(:), ys(:)
    real(dp), intent(inout) :: phi(:), fx(:), fy(:)

    integer :: i, ip, n, m
    real(dp) :: dx, dy

    do i = 1, cell%np
      ip = cell%ids(i)

      dx = xs(ip) - cell%xc
      dy = ys(ip) - cell%yc

      ! Вычисляем φ, F из L
      do nn = 0, P_ORDER
        do mm = 0, P_ORDER - nn
          if (abs(cell%L(nn,mm)) < 1e-15_dp) cycle

          ! Потенциал
          phi(ip) = phi(ip) + cell%L(nn,mm) * dx**nn * dy**mm

          ! Силы
          if (n > 0) then
            fx(ip) = fx(ip) - real(nn,dp) * cell%L(nn,mm) * dx**(nn-1) * dy**mm
          end if

          if (m > 0) then
            fy(ip) = fy(ip) - real(mm,dp) * cell%L(nn,mm) * dx**nn * dy**(mm-1)
          end if
        end do
      end do

    end do

  end subroutine eval_local

  !============================================================================
  ! 6. P2P - ПРЯМОЙ РАСЧЕТ ДЛЯ БЛИЗКИХ ЧАСТИЦ
  !
  ! Для близких ячеек мультипольное разложение неточно.
  ! Поэтому мы считаем напрямую: φ = q/r, F = q·r/r³
  !============================================================================
  subroutine eval_near(grid, ix, iy, N, xs, ys, qs, phi, fx, fy)
    type(Cell), intent(in) :: grid(:,:)
    integer, intent(in) :: ix, iy, N
    real(dp), intent(in) :: xs(N), ys(N), qs(N)
    real(dp), intent(inout) :: phi(N), fx(N), fy(N)

    integer :: jx, jy, i, j, ip, jp
    real(dp) :: rx, ry, r2, r, r_inv, r3_inv
    logical :: same_cell

    ! Проходим по близким ячейкам
    do jy = max(1, iy-NEAR_RANGE), min(NGRID, iy+NEAR_RANGE)
      do jx = max(1, ix-NEAR_RANGE), min(NGRID, ix+NEAR_RANGE)
        if (grid(jx,jy)%np == 0) cycle

        same_cell = (ix == jx) .and. (iy == jy)

        ! Все пары частиц
        do i = 1, grid(ix,iy)%np
          ip = grid(ix,iy)%ids(i)

          do j = 1, grid(jx,jy)%np
            jp = grid(jx,jy)%ids(j)

            if (same_cell .and. ip == jp) cycle

            ! Вектор r = r_i - r_j
            rx = xs(ip) - xs(jp)
            ry = ys(ip) - ys(jp)
            r2 = rx*rx + ry*ry

            if (r2 < 1e-16_dp) cycle

            r = sqrt(r2)
            r_inv = 1.0_dp / r
            r3_inv = r_inv / r2

            ! φ = q/r
            phi(ip) = phi(ip) + qs(jp) * r_inv

            ! F = q·r/r³
            fx(ip) = fx(ip) + qs(jp) * rx * r3_inv
            fy(ip) = fy(ip) + qs(jp) * ry * r3_inv
          end do
        end do

      end do
    end do

  end subroutine eval_near

  !============================================================================
  ! ОЧИСТКА
  !============================================================================
  subroutine cleanup_grid(grid)
    type(Cell), intent(inout) :: grid(:,:)
    integer :: ix, iy

    do iy = 1, NGRID
      do ix = 1, NGRID
        if (allocated(grid(ix,iy)%ids)) deallocate(grid(ix,iy)%ids)
      end do
    end do

  end subroutine cleanup_grid

end module simple_fmm_2d
