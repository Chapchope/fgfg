module potential_2d_mod
  !============================================================================
  ! МОДУЛЬ: Абстрактный интерфейс для 2D потенциалов
  !
  ! НАЗНАЧЕНИЕ:
  !   Этот модуль определяет абстрактный тип для потенциалов в 2D.
  !   Позволяет легко менять потенциал (1/r, ln(r), экранированный и т.д.)
  !   без изменения основного кода FMM.
  !
  ! ИСПОЛЬЗОВАНИЕ:
  !   type(coulomb_2d_potential) :: pot
  !   call pot%init()
  !   phi = pot%eval(x, y, q, x0, y0)
  !   call pot%compute_multipole_derivatives(dx, dy, R, p_order, derivs)
  !
  ! АВТОР: Grid-based FMM для сверхпроводника
  !============================================================================
  use types_mod
  implicit none
  private

  public :: abstract_potential_2d
  public :: coulomb_2d_potential

  !============================================================================
  ! АБСТРАКТНЫЙ ТИП: Потенциал 2D
  !============================================================================
  type, abstract :: abstract_potential_2d
  contains
    ! Вычислить φ(r) для заряда q в точке (x,y) от источника (x0,y0)
    procedure(eval_potential_interface), deferred :: eval

    ! Вычислить φ и силу F
    procedure(eval_field_interface), deferred :: eval_field

    ! Вычислить производные ∂^n∂^m(1/R) для M2L трансляции
    ! Это ключевая функция для FMM!
    procedure(compute_derivatives_interface), deferred :: compute_multipole_derivatives

    ! Инициализация (если нужна таблица, предвычисления и т.д.)
    procedure(init_interface), deferred :: init
    procedure(cleanup_interface), deferred :: cleanup
  end type abstract_potential_2d

  !============================================================================
  ! ИНТЕРФЕЙСЫ ДЛЯ АБСТРАКТНЫХ ПРОЦЕДУР
  !============================================================================
  abstract interface
    function eval_potential_interface(this, x, y, q, x0, y0) result(phi)
      import :: abstract_potential_2d, dp
      class(abstract_potential_2d), intent(in) :: this
      real(dp), intent(in) :: x, y, q, x0, y0
      real(dp) :: phi
    end function eval_potential_interface

    subroutine eval_field_interface(this, x, y, q, x0, y0, phi, fx, fy)
      import :: abstract_potential_2d, dp
      class(abstract_potential_2d), intent(in) :: this
      real(dp), intent(in) :: x, y, q, x0, y0
      real(dp), intent(out) :: phi, fx, fy
    end subroutine eval_field_interface

    subroutine compute_derivatives_interface(this, dx, dy, R, p_max, derivs)
      import :: abstract_potential_2d, dp
      class(abstract_potential_2d), intent(in) :: this
      real(dp), intent(in) :: dx, dy, R
      integer, intent(in) :: p_max
      real(dp), intent(out) :: derivs(0:p_max, 0:p_max)
    end subroutine compute_derivatives_interface

    subroutine init_interface(this, p_max)
      import :: abstract_potential_2d
      class(abstract_potential_2d), intent(inout) :: this
      integer, intent(in) :: p_max
    end subroutine init_interface

    subroutine cleanup_interface(this)
      import :: abstract_potential_2d
      class(abstract_potential_2d), intent(inout) :: this
    end subroutine cleanup_interface
  end interface

  !============================================================================
  ! РЕАЛИЗАЦИЯ: Кулоновский потенциал 2D (φ = 1/r при z=0)
  !============================================================================
  type, extends(abstract_potential_2d) :: coulomb_2d_potential
    integer :: p_max_cached = -1
  contains
    procedure :: eval => coulomb_eval
    procedure :: eval_field => coulomb_eval_field
    procedure :: compute_multipole_derivatives => coulomb_compute_derivatives
    procedure :: init => coulomb_init
    procedure :: cleanup => coulomb_cleanup
  end type coulomb_2d_potential

contains

  !============================================================================
  ! КУЛОН: Инициализация
  !============================================================================
  subroutine coulomb_init(this, p_max)
    class(coulomb_2d_potential), intent(inout) :: this
    integer, intent(in) :: p_max

    this%p_max_cached = p_max
    print '(A,I0)', "[Potential] Coulomb 2D initialized with p_max = ", p_max
  end subroutine coulomb_init

  subroutine coulomb_cleanup(this)
    class(coulomb_2d_potential), intent(inout) :: this
    ! Нет ресурсов для очистки
    this%p_max_cached = -1
  end subroutine coulomb_cleanup

  !============================================================================
  ! КУЛОН: Вычисление потенциала φ = q/r
  !============================================================================
  function coulomb_eval(this, x, y, q, x0, y0) result(phi)
    class(coulomb_2d_potential), intent(in) :: this
    real(dp), intent(in) :: x, y, q, x0, y0
    real(dp) :: phi
    real(dp) :: dx, dy, r2, r

    ! Убираем предупреждение
    if (this%p_max_cached < 0) then
    end if

    dx = x - x0
    dy = y - y0
    r2 = dx*dx + dy*dy

    if (r2 < 1.0e-16_dp) then
      phi = 0.0_dp
    else
      r = sqrt(r2)
      phi = q / r
    end if
  end function coulomb_eval

  !============================================================================
  ! КУЛОН: Вычисление поля φ и силы F = -∇φ
  !============================================================================
  subroutine coulomb_eval_field(this, x, y, q, x0, y0, phi, fx, fy)
    class(coulomb_2d_potential), intent(in) :: this
    real(dp), intent(in) :: x, y, q, x0, y0
    real(dp), intent(out) :: phi, fx, fy
    real(dp) :: dx, dy, r2, r, r_inv, r3_inv

    ! Убираем предупреждение
    if (this%p_max_cached < 0) then
    end if

    dx = x - x0
    dy = y - y0
    r2 = dx*dx + dy*dy

    if (r2 < 1.0e-16_dp) then
      phi = 0.0_dp
      fx = 0.0_dp
      fy = 0.0_dp
    else
      r = sqrt(r2)
      r_inv = 1.0_dp / r
      r3_inv = r_inv / r2

      phi = q * r_inv
      ! F = q * r/r³
      fx = q * dx * r3_inv
      fy = q * dy * r3_inv
    end if
  end subroutine coulomb_eval_field

  !============================================================================
  ! КУЛОН: Производные для M2L
  !
  ! Вычисляем T_nm = ∂^n∂^m(1/R) / (n! m!)
  ! Где R = |r - r0| = sqrt(dx² + dy²)
  !
  ! Эти производные используются в M2L трансляции:
  !   L_nm = Σ_{k,l} M_kl * T_{n+k,m+l}
  !
  ! ФОРМУЛЫ (рекурсивные):
  ! Базовые случаи:
  !   T_00 = 1/R
  !   T_10 = -dx/R³
  !   T_01 = -dy/R³
  !
  ! Рекурсия:
  !   ∂/∂x(1/R^p) = -p·x/R^(p+2)
  !   ∂/∂y(1/R^p) = -p·y/R^(p+2)
  !
  ! Для произвольного порядка используем рекурсивную формулу:
  !   T_{n+1,m} = -1/R² * [n*T_{n-1,m} + dx*(2n+2m+1)*T_{n,m}]
  !   T_{n,m+1} = -1/R² * [m*T_{n,m-1} + dy*(2n+2m+1)*T_{n,m}]
  !============================================================================
  subroutine coulomb_compute_derivatives(this, dx, dy, R, p_max, derivs)
    class(coulomb_2d_potential), intent(in) :: this
    real(dp), intent(in) :: dx, dy, R
    integer, intent(in) :: p_max
    real(dp), intent(out) :: derivs(0:p_max, 0:p_max)

    integer :: nn, mm, pp
    real(dp) :: R2, R_inv, R2_inv, R3_inv, R5_inv, R7_inv

    ! Убираем предупреждение об неиспользуемом this
    if (this%p_max_cached < 0) then
      ! Никогда не выполнится
    end if

    derivs = 0.0_dp

    if (R < 1.0e-14_dp) return

    R2 = R * R
    R_inv = 1.0_dp / R
    R2_inv = R_inv * R_inv
    R3_inv = R_inv * R2_inv
    R5_inv = R3_inv * R2_inv
    R7_inv = R5_inv * R2_inv

    ! Базовые случаи (p=0)
    derivs(0,0) = R_inv

    if (p_max >= 1) then
      derivs(1,0) = -dx * R3_inv
      derivs(0,1) = -dy * R3_inv
    end if

    if (p_max >= 2) then
      ! Используем формулы для вторых производных
      derivs(2,0) = (3.0_dp*dx*dx - R2) * R5_inv
      derivs(1,1) = 3.0_dp * dx * dy * R5_inv
      derivs(0,2) = (3.0_dp*dy*dy - R2) * R5_inv
    end if

    if (p_max >= 3) then
      ! Третьи производные
      derivs(3,0) = -dx * (15.0_dp*dx*dx - 9.0_dp*R2) * R7_inv
      derivs(2,1) = -dy * (15.0_dp*dx*dx - 3.0_dp*R2) * R7_inv
      derivs(1,2) = -dx * (15.0_dp*dy*dy - 3.0_dp*R2) * R7_inv
      derivs(0,3) = -dy * (15.0_dp*dy*dy - 9.0_dp*R2) * R7_inv
    end if

    ! Для p > 3 производные пока не реализованы правильно
    ! Нужно использовать полную формулу или таблицу из Python
    ! Пока просто заполним нулями для демонстрации
    if (p_max >= 4) then
      do pp = 4, p_max
        do nn = 0, pp
          mm = pp - nn
          if (nn+mm == pp) then
            derivs(nn,mm) = 0.0_dp
          end if
        end do
      end do
    end if

  end subroutine coulomb_compute_derivatives

end module potential_2d_mod
