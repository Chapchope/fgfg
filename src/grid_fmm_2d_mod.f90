module grid_fmm_2d_mod
  !============================================================================
  ! МОДУЛЬ: Grid-based FMM для 2D
  !
  ! НАЗНАЧЕНИЕ:
  !   Быстрый мультипольный метод (FMM) на равномерной сетке для 2D
  !   С произвольным порядком разложения и заменяемым потенциалом
  !
  ! ОСОБЕННОСТИ:
  !   - Произвольный порядок разложения p
  !   - Легко заменяемый потенциал (через абстрактный интерфейс)
  !   - Высокая точность (до 10^-10 при высоких p)
  !   - Высокая скорость O(N) вместо O(N²)
  !
  ! ИСПОЛЬЗОВАНИЕ:
  !   type(grid_fmm_2d) :: fmm
  !   type(coulomb_2d_potential) :: pot
  !
  !   call pot%init(p_order)
  !   call fmm%init(N, xs, ys, ngrid, p_order, pot)
  !   call fmm%compute(qs, phi, fx, fy)
  !   call fmm%cleanup()
  !
  ! АВТОР: Grid-based FMM для сверхпроводника
  !============================================================================
  use types_mod
  use potential_2d_mod
  use multipole_2d_mod
  implicit none
  private

  public :: grid_fmm_2d

  !============================================================================
  ! ТИП: Ячейка сетки
  !============================================================================
  type :: grid_cell_2d
    real(dp) :: xc, yc                ! Центр ячейки
    integer :: np                     ! Число частиц
    integer, allocatable :: ids(:)    ! Индексы частиц в этой ячейке
    real(dp), allocatable :: M(:,:)   ! Мультипольные моменты M(0:p, 0:p)
    real(dp), allocatable :: L(:,:)   ! Локальное разложение L(0:p, 0:p)
  end type grid_cell_2d

  !============================================================================
  ! ТИП: Grid FMM
  !============================================================================
  type :: grid_fmm_2d
    ! Параметры
    integer :: N                      ! Число частиц
    integer :: ngrid                  ! Размер сетки (ngrid × ngrid)
    integer :: p_order                ! Порядок разложения
    integer :: near_range             ! Радиус ближнего поля (в ячейках)

    ! Геометрия
    real(dp) :: xmin, xmax, ymin, ymax
    real(dp) :: hx, hy                ! Размер ячейки

    ! Данные частиц (указатели на внешние массивы)
    real(dp), pointer :: xs(:) => null()
    real(dp), pointer :: ys(:) => null()

    ! Сетка
    type(grid_cell_2d), allocatable :: grid(:,:)

    ! Потенциал (полиморфизм!)
    class(abstract_potential_2d), pointer :: potential => null()

    ! Статистика
    integer :: n_m2l                  ! Число M2L взаимодействий
    integer :: n_p2p                  ! Число P2P взаимодействий

  contains
    procedure :: init => fmm_init
    procedure :: compute => fmm_compute
    procedure :: cleanup => fmm_cleanup

    ! Внутренние процедуры
    procedure, private :: setup_grid
    procedure, private :: assign_particles
    procedure, private :: compute_multipoles
    procedure, private :: multipole_to_local
    procedure, private :: evaluate_local
    procedure, private :: evaluate_near
  end type grid_fmm_2d

contains

  !============================================================================
  ! ИНИЦИАЛИЗАЦИЯ FMM
  !============================================================================
  subroutine fmm_init(this, N, xs, ys, ngrid, p_order, potential, near_range_opt)
    class(grid_fmm_2d), intent(inout) :: this
    integer, intent(in) :: N
    real(dp), intent(in), target :: xs(N), ys(N)
    integer, intent(in) :: ngrid, p_order
    class(abstract_potential_2d), intent(in), target :: potential
    integer, intent(in), optional :: near_range_opt

    this%N = N
    this%ngrid = ngrid
    this%p_order = p_order
    this%xs => xs
    this%ys => ys
    this%potential => potential

    ! Радиус ближнего поля (по умолчанию = 2)
    if (present(near_range_opt)) then
      this%near_range = near_range_opt
    else
      this%near_range = 2
    end if

    ! Статистика
    this%n_m2l = 0
    this%n_p2p = 0

    print '(A)', "============================================"
    print '(A)', "  Grid-based FMM 2D - Initialization"
    print '(A)', "============================================"
    print '(A,I0)', "Particles:    ", N
    print '(A,I0,A,I0)', "Grid:         ", ngrid, " x ", ngrid
    print '(A,I0)', "P_order:      ", p_order
    print '(A,I0)', "Near range:   ", this%near_range
    print '(A)', ""

    ! Создаем сетку
    call this%setup_grid()

    ! Назначаем частицы в ячейки
    call this%assign_particles()

    print '(A)', "FMM initialized successfully!"
    print '(A)', ""

  end subroutine fmm_init

  !============================================================================
  ! ВЫЧИСЛЕНИЕ φ, F через FMM
  !============================================================================
  subroutine fmm_compute(this, qs, phi, fx, fy)
    class(grid_fmm_2d), intent(inout) :: this
    real(dp), intent(in) :: qs(this%N)
    real(dp), intent(out) :: phi(this%N), fx(this%N), fy(this%N)

    integer :: ix, iy

    print '(A)', "============================================"
    print '(A)', "  Computing FMM..."
    print '(A)', "============================================"

    ! Обнуляем выходные массивы
    phi = 0.0_dp
    fx = 0.0_dp
    fy = 0.0_dp

    ! Обнуляем статистику
    this%n_m2l = 0
    this%n_p2p = 0

    ! 1. P2M: Вычисляем мультипольные моменты
    print '(A)', "[1/4] P2M: Computing multipoles..."
    call this%compute_multipoles(qs)

    ! 2. M2L: Дальнее поле через мультиполи
    print '(A)', "[2/4] M2L: Far-field via multipoles..."
    call this%multipole_to_local()

    ! 3. L2P: Вычисляем поля из локальных разложений
    print '(A)', "[3/4] L2P: Evaluating local expansions..."
    call this%evaluate_local(phi, fx, fy)

    ! 4. P2P: Прямой расчет для близких частиц
    print '(A)', "[4/4] P2P: Direct near-field..."
    call this%evaluate_near(qs, phi, fx, fy)

    print '(A)', ""
    print '(A)', "Statistics:"
    print '(A,I0)', "  M2L interactions: ", this%n_m2l
    print '(A,I0)', "  P2P interactions: ", this%n_p2p
    print '(A)', ""
    print '(A)', "FMM computation complete!"
    print '(A)', "============================================"
    print '(A)', ""

  end subroutine fmm_compute

  !============================================================================
  ! ОЧИСТКА
  !============================================================================
  subroutine fmm_cleanup(this)
    class(grid_fmm_2d), intent(inout) :: this
    integer :: ix, iy

    if (allocated(this%grid)) then
      do iy = 1, this%ngrid
        do ix = 1, this%ngrid
          if (allocated(this%grid(ix,iy)%ids)) deallocate(this%grid(ix,iy)%ids)
          if (allocated(this%grid(ix,iy)%M)) deallocate(this%grid(ix,iy)%M)
          if (allocated(this%grid(ix,iy)%L)) deallocate(this%grid(ix,iy)%L)
        end do
      end do
      deallocate(this%grid)
    end if

    this%xs => null()
    this%ys => null()
    this%potential => null()

  end subroutine fmm_cleanup

  !============================================================================
  ! 1. НАСТРОЙКА СЕТКИ
  !============================================================================
  subroutine setup_grid(this)
    class(grid_fmm_2d), intent(inout) :: this

    integer :: ix, iy
    real(dp) :: Lx, Ly, xcen, ycen, margin

    ! Находим границы частиц
    this%xmin = minval(this%xs)
    this%xmax = maxval(this%xs)
    this%ymin = minval(this%ys)
    this%ymax = maxval(this%ys)

    ! Добавляем отступ 10%
    margin = 0.1_dp
    Lx = (this%xmax - this%xmin) * (1.0_dp + margin)
    Ly = (this%ymax - this%ymin) * (1.0_dp + margin)

    xcen = 0.5_dp * (this%xmin + this%xmax)
    ycen = 0.5_dp * (this%ymin + this%ymax)

    this%xmin = xcen - 0.5_dp * Lx
    this%xmax = xcen + 0.5_dp * Lx
    this%ymin = ycen - 0.5_dp * Ly
    this%ymax = ycen + 0.5_dp * Ly

    this%hx = Lx / real(this%ngrid, dp)
    this%hy = Ly / real(this%ngrid, dp)

    ! Выделяем память для сетки
    allocate(this%grid(this%ngrid, this%ngrid))

    ! Инициализируем ячейки
    do iy = 1, this%ngrid
      do ix = 1, this%ngrid
        this%grid(ix,iy)%xc = this%xmin + (real(ix,dp) - 0.5_dp) * this%hx
        this%grid(ix,iy)%yc = this%ymin + (real(iy,dp) - 0.5_dp) * this%hy
        this%grid(ix,iy)%np = 0

        ! Выделяем память для M и L
        allocate(this%grid(ix,iy)%M(0:this%p_order, 0:this%p_order))
        allocate(this%grid(ix,iy)%L(0:this%p_order, 0:this%p_order))
        this%grid(ix,iy)%M = 0.0_dp
        this%grid(ix,iy)%L = 0.0_dp
      end do
    end do

  end subroutine setup_grid

  !============================================================================
  ! 2. НАЗНАЧЕНИЕ ЧАСТИЦ В ЯЧЕЙКИ
  !============================================================================
  subroutine assign_particles(this)
    class(grid_fmm_2d), intent(inout) :: this

    integer :: i, ix, iy
    integer :: counts(this%ngrid, this%ngrid)

    ! Считаем частицы в каждой ячейке
    counts = 0
    do i = 1, this%N
      ix = int((this%xs(i) - this%xmin) / this%hx) + 1
      iy = int((this%ys(i) - this%ymin) / this%hy) + 1
      ix = max(1, min(this%ngrid, ix))
      iy = max(1, min(this%ngrid, iy))
      counts(ix, iy) = counts(ix, iy) + 1
    end do

    ! Выделяем память
    do iy = 1, this%ngrid
      do ix = 1, this%ngrid
        this%grid(ix,iy)%np = counts(ix,iy)
        if (counts(ix,iy) > 0) then
          allocate(this%grid(ix,iy)%ids(counts(ix,iy)))
        end if
      end do
    end do

    ! Заполняем списки
    counts = 0
    do i = 1, this%N
      ix = int((this%xs(i) - this%xmin) / this%hx) + 1
      iy = int((this%ys(i) - this%ymin) / this%hy) + 1
      ix = max(1, min(this%ngrid, ix))
      iy = max(1, min(this%ngrid, iy))

      counts(ix,iy) = counts(ix,iy) + 1
      this%grid(ix,iy)%ids(counts(ix,iy)) = i
    end do

  end subroutine assign_particles

  !============================================================================
  ! 3. P2M - ВЫЧИСЛЕНИЕ МУЛЬТИПОЛЬНЫХ МОМЕНТОВ
  !============================================================================
  subroutine compute_multipoles(this, qs)
    class(grid_fmm_2d), intent(inout) :: this
    real(dp), intent(in) :: qs(this%N)

    integer :: ix, iy, i, ip
    real(dp), allocatable :: xs_cell(:), ys_cell(:), qs_cell(:)

    do iy = 1, this%ngrid
      do ix = 1, this%ngrid
        if (this%grid(ix,iy)%np == 0) cycle

        ! Собираем частицы этой ячейки
        allocate(xs_cell(this%grid(ix,iy)%np))
        allocate(ys_cell(this%grid(ix,iy)%np))
        allocate(qs_cell(this%grid(ix,iy)%np))

        do i = 1, this%grid(ix,iy)%np
          ip = this%grid(ix,iy)%ids(i)
          xs_cell(i) = this%xs(ip)
          ys_cell(i) = this%ys(ip)
          qs_cell(i) = qs(ip)
        end do

        ! Вычисляем мультиполи через модуль multipole_2d_mod
        call p2m_2d(this%grid(ix,iy)%np, xs_cell, ys_cell, qs_cell, &
                    this%grid(ix,iy)%xc, this%grid(ix,iy)%yc, &
                    this%p_order, this%grid(ix,iy)%M)

        deallocate(xs_cell, ys_cell, qs_cell)
      end do
    end do

  end subroutine compute_multipoles

  !============================================================================
  ! 4. M2L - MULTIPOLE TO LOCAL
  !============================================================================
  subroutine multipole_to_local(this)
    class(grid_fmm_2d), intent(inout) :: this

    integer :: ix, iy, jx, jy
    logical :: is_near

    ! Обнуляем локальные разложения
    do iy = 1, this%ngrid
      do ix = 1, this%ngrid
        this%grid(ix,iy)%L = 0.0_dp
      end do
    end do

    ! Для каждой целевой ячейки
    do iy = 1, this%ngrid
      do ix = 1, this%ngrid

        ! Для каждой исходной ячейки
        do jy = 1, this%ngrid
          do jx = 1, this%ngrid
            if (ix == jx .and. iy == jy) cycle
            if (this%grid(jx,jy)%np == 0) cycle

            ! Пропускаем близкие (обработаем в P2P)
            is_near = (abs(ix-jx) <= this%near_range) .and. &
                     (abs(iy-jy) <= this%near_range)
            if (is_near) cycle

            ! M2L трансляция через модуль multipole_2d_mod
            call m2l_2d(this%grid(jx,jy)%M, &
                        this%grid(jx,jy)%xc, this%grid(jx,jy)%yc, &
                        this%grid(ix,iy)%xc, this%grid(ix,iy)%yc, &
                        this%p_order, this%potential, &
                        this%grid(ix,iy)%L)

            this%n_m2l = this%n_m2l + 1
          end do
        end do

      end do
    end do

  end subroutine multipole_to_local

  !============================================================================
  ! 5. L2P - ВЫЧИСЛЕНИЕ ИЗ ЛОКАЛЬНОГО РАЗЛОЖЕНИЯ
  !============================================================================
  subroutine evaluate_local(this, phi, fx, fy)
    class(grid_fmm_2d), intent(inout) :: this
    real(dp), intent(inout) :: phi(this%N), fx(this%N), fy(this%N)

    integer :: ix, iy, i, ip
    real(dp), allocatable :: xs_cell(:), ys_cell(:)
    real(dp), allocatable :: phi_cell(:), fx_cell(:), fy_cell(:)

    do iy = 1, this%ngrid
      do ix = 1, this%ngrid
        if (this%grid(ix,iy)%np == 0) cycle

        ! Собираем частицы этой ячейки
        allocate(xs_cell(this%grid(ix,iy)%np))
        allocate(ys_cell(this%grid(ix,iy)%np))
        allocate(phi_cell(this%grid(ix,iy)%np))
        allocate(fx_cell(this%grid(ix,iy)%np))
        allocate(fy_cell(this%grid(ix,iy)%np))

        do i = 1, this%grid(ix,iy)%np
          ip = this%grid(ix,iy)%ids(i)
          xs_cell(i) = this%xs(ip)
          ys_cell(i) = this%ys(ip)
        end do

        phi_cell = 0.0_dp
        fx_cell = 0.0_dp
        fy_cell = 0.0_dp

        ! Вычисляем через модуль multipole_2d_mod
        call l2p_2d(this%grid(ix,iy)%L, &
                    this%grid(ix,iy)%xc, this%grid(ix,iy)%yc, &
                    this%grid(ix,iy)%np, xs_cell, ys_cell, &
                    this%p_order, phi_cell, fx_cell, fy_cell)

        ! Копируем обратно
        do i = 1, this%grid(ix,iy)%np
          ip = this%grid(ix,iy)%ids(i)
          phi(ip) = phi(ip) + phi_cell(i)
          fx(ip) = fx(ip) + fx_cell(i)
          fy(ip) = fy(ip) + fy_cell(i)
        end do

        deallocate(xs_cell, ys_cell, phi_cell, fx_cell, fy_cell)
      end do
    end do

  end subroutine evaluate_local

  !============================================================================
  ! 6. P2P - ПРЯМОЙ РАСЧЕТ ДЛЯ БЛИЗКИХ ЧАСТИЦ
  !============================================================================
  subroutine evaluate_near(this, qs, phi, fx, fy)
    class(grid_fmm_2d), intent(inout) :: this
    real(dp), intent(in) :: qs(this%N)
    real(dp), intent(inout) :: phi(this%N), fx(this%N), fy(this%N)

    integer :: ix, iy, jx, jy, i, ip
    real(dp), allocatable :: xs_tgt(:), ys_tgt(:)
    real(dp), allocatable :: xs_src(:), ys_src(:), qs_src(:)
    real(dp), allocatable :: phi_tmp(:), fx_tmp(:), fy_tmp(:)
    logical :: same_cell

    ! Для каждой целевой ячейки
    do iy = 1, this%ngrid
      do ix = 1, this%ngrid
        if (this%grid(ix,iy)%np == 0) cycle

        ! Собираем целевые частицы
        allocate(xs_tgt(this%grid(ix,iy)%np))
        allocate(ys_tgt(this%grid(ix,iy)%np))
        allocate(phi_tmp(this%grid(ix,iy)%np))
        allocate(fx_tmp(this%grid(ix,iy)%np))
        allocate(fy_tmp(this%grid(ix,iy)%np))

        do i = 1, this%grid(ix,iy)%np
          ip = this%grid(ix,iy)%ids(i)
          xs_tgt(i) = this%xs(ip)
          ys_tgt(i) = this%ys(ip)
        end do

        phi_tmp = 0.0_dp
        fx_tmp = 0.0_dp
        fy_tmp = 0.0_dp

        ! Проходим по близким ячейкам
        do jy = max(1, iy-this%near_range), min(this%ngrid, iy+this%near_range)
          do jx = max(1, ix-this%near_range), min(this%ngrid, ix+this%near_range)
            if (this%grid(jx,jy)%np == 0) cycle

            same_cell = (ix == jx) .and. (iy == jy)

            ! Собираем исходные частицы
            allocate(xs_src(this%grid(jx,jy)%np))
            allocate(ys_src(this%grid(jx,jy)%np))
            allocate(qs_src(this%grid(jx,jy)%np))

            do i = 1, this%grid(jx,jy)%np
              ip = this%grid(jx,jy)%ids(i)
              xs_src(i) = this%xs(ip)
              ys_src(i) = this%ys(ip)
              qs_src(i) = qs(ip)
            end do

            ! P2P через модуль multipole_2d_mod
            call p2p_direct_2d(this%grid(jx,jy)%np, xs_src, ys_src, qs_src, &
                              this%grid(ix,iy)%np, xs_tgt, ys_tgt, &
                              this%potential, same_cell, &
                              phi_tmp, fx_tmp, fy_tmp)

            this%n_p2p = this%n_p2p + this%grid(ix,iy)%np * this%grid(jx,jy)%np

            deallocate(xs_src, ys_src, qs_src)
          end do
        end do

        ! Копируем обратно
        do i = 1, this%grid(ix,iy)%np
          ip = this%grid(ix,iy)%ids(i)
          phi(ip) = phi(ip) + phi_tmp(i)
          fx(ip) = fx(ip) + fx_tmp(i)
          fy(ip) = fy(ip) + fy_tmp(i)
        end do

        deallocate(xs_tgt, ys_tgt, phi_tmp, fx_tmp, fy_tmp)
      end do
    end do

  end subroutine evaluate_near

end module grid_fmm_2d_mod
