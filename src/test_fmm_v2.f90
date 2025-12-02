program test_fmm_v2
  !============================================================================
  ! ТЕСТ: Grid-based FMM 2D с произвольным порядком и заменяемым потенциалом
  !
  ! Проверяем:
  !   1. Точность при разных p_order (2, 4, 6, 8, 10)
  !   2. Сравнение с прямым расчетом
  !   3. Скорость
  !
  ! АВТОР: Grid-based FMM для сверхпроводника
  !============================================================================
  use types_mod
  use potential_2d_mod
  use grid_fmm_2d_mod
  implicit none

  integer, parameter :: N = 1000                ! Число частиц
  integer, parameter :: NGRID = 30              ! Размер сетки
  integer, parameter :: P_ORDER = 3             ! Порядок разложения (МЕНЯЙТЕ!)
  integer, parameter :: NEAR_RANGE = 5          ! Радиус ближнего поля

  real(dp), allocatable :: xs(:), ys(:), qs(:)
  real(dp), allocatable :: phi_fmm(:), fx_fmm(:), fy_fmm(:)
  real(dp), allocatable :: phi_direct(:), fx_direct(:), fy_direct(:)

  type(grid_fmm_2d) :: fmm
  type(coulomb_2d_potential), target :: pot

  real(dp) :: t1, t2, time_fmm, time_direct
  real(dp) :: err_phi, err_fx, err_fy
  integer :: i

  print '(A)', ""
  print '(A)', "============================================"
  print '(A)', "  TEST: Grid-based FMM 2D"
  print '(A)', "============================================"
  print '(A)', ""

  ! Выделяем память
  allocate(xs(N), ys(N), qs(N))
  allocate(phi_fmm(N), fx_fmm(N), fy_fmm(N))
  allocate(phi_direct(N), fx_direct(N), fy_direct(N))

  ! Генерируем случайные частицы
  call random_seed()
  call random_number(xs)
  call random_number(ys)
  call random_number(qs)

  ! Нормируем
  xs = xs * 10.0_dp - 5.0_dp  ! [-5, 5]
  ys = ys * 10.0_dp - 5.0_dp
  qs = (qs - 0.5_dp) * 2.0_dp ! [-1, 1]

  print '(A,I0)', "Generated ", N, " random particles"
  print '(A)', ""

  !============================================================================
  ! 1. ПРЯМОЙ РАСЧЕТ (для проверки точности)
  !============================================================================
  print '(A)', "--------------------------------------------"
  print '(A)', " Computing direct (exact) solution..."
  print '(A)', "--------------------------------------------"

  call cpu_time(t1)
  call compute_direct(N, xs, ys, qs, pot, phi_direct, fx_direct, fy_direct)
  call cpu_time(t2)
  time_direct = t2 - t1

  print '(A,F10.6,A)', "Direct computation time: ", time_direct, " s"
  print '(A,ES12.4)', "Sample phi_direct(1):    ", phi_direct(1)
  print '(A)', ""

  !============================================================================
  ! 2. FMM РАСЧЕТ
  !============================================================================
  print '(A)', "--------------------------------------------"
  print '(A)', " Computing FMM solution..."
  print '(A)', "--------------------------------------------"

  ! Инициализируем потенциал
  call pot%init(P_ORDER)

  ! Инициализируем FMM
  call fmm%init(N, xs, ys, NGRID, P_ORDER, pot, NEAR_RANGE)

  ! Вычисляем через FMM
  call cpu_time(t1)
  call fmm%compute(qs, phi_fmm, fx_fmm, fy_fmm)
  call cpu_time(t2)
  time_fmm = t2 - t1

  print '(A,F10.6,A)', "FMM computation time:    ", time_fmm, " s"
  print '(A,ES12.4)', "Sample phi_fmm(1):       ", phi_fmm(1)
  print '(A)', ""

  ! Очищаем FMM
  call fmm%cleanup()
  call pot%cleanup()

  !============================================================================
  ! 3. ПРОВЕРКА ТОЧНОСТИ
  !============================================================================
  print '(A)', "============================================"
  print '(A)', "  Accuracy Check"
  print '(A)', "============================================"

  err_phi = 0.0_dp
  err_fx = 0.0_dp
  err_fy = 0.0_dp

  do i = 1, N
    if (abs(phi_direct(i)) > 1.0e-10_dp) then
      err_phi = max(err_phi, abs(phi_fmm(i) - phi_direct(i)) / abs(phi_direct(i)))
    end if

    if (abs(fx_direct(i)) > 1.0e-10_dp) then
      err_fx = max(err_fx, abs(fx_fmm(i) - fx_direct(i)) / abs(fx_direct(i)))
    end if

    if (abs(fy_direct(i)) > 1.0e-10_dp) then
      err_fy = max(err_fy, abs(fy_fmm(i) - fy_direct(i)) / abs(fy_direct(i)))
    end if
  end do

  print '(A,ES12.4)', "Max relative error (phi): ", err_phi
  print '(A,ES12.4)', "Max relative error (Fx):  ", err_fx
  print '(A,ES12.4)', "Max relative error (Fy):  ", err_fy
  print '(A)', ""

  ! Оценка успешности
  print '(A)', "--------------------------------------------"
  if (err_phi < 1.0e-6_dp) then
    print '(A)', "✓ Excellent! Error < 10^-6"
  else if (err_phi < 1.0e-4_dp) then
    print '(A)', "✓ Good! Error < 10^-4"
  else if (err_phi < 1.0e-2_dp) then
    print '(A)', "○ Acceptable. Error < 1%"
  else
    print '(A)', "✗ Poor accuracy. Need higher p_order or larger near_range"
  end if
  print '(A)', "--------------------------------------------"
  print '(A)', ""

  !============================================================================
  ! 4. СРАВНЕНИЕ СКОРОСТИ
  !============================================================================
  print '(A)', "============================================"
  print '(A)', "  Performance"
  print '(A)', "============================================"
  print '(A,F10.6,A)', "Direct time:  ", time_direct, " s"
  print '(A,F10.6,A)', "FMM time:     ", time_fmm, " s"
  print '(A,F10.2,A)', "Speedup:      ", time_direct / time_fmm, "×"
  print '(A)', ""

  if (time_fmm < time_direct) then
    print '(A)', "✓ FMM is faster!"
  else
    print '(A)', "○ For N=1000, direct may be faster."
    print '(A)', "  Try N=10000 to see FMM advantage!"
  end if
  print '(A)', ""

  !============================================================================
  ! 5. РЕКОМЕНДАЦИИ
  !============================================================================
  print '(A)', "============================================"
  print '(A)', "  Recommendations"
  print '(A)', "============================================"

  if (err_phi > 1.0e-6_dp) then
    print '(A)', "To improve accuracy:"
    print '(A)', "  1. Increase P_ORDER (current: ", P_ORDER, ")"
    print '(A)', "  2. Increase NEAR_RANGE (current: ", NEAR_RANGE, ")"
    print '(A)', "  3. Check potential derivatives implementation"
    print '(A)', ""
  else
    print '(A)', "Accuracy is excellent!"
    print '(A)', ""
  end if

  if (time_fmm >= time_direct) then
    print '(A)', "To see FMM speed advantage:"
    print '(A)', "  - Increase N to 10,000 or 100,000"
    print '(A)', "  - FMM scales as O(N), direct as O(N²)"
    print '(A)', ""
  end if

  print '(A)', "============================================"
  print '(A)', ""

  deallocate(xs, ys, qs)
  deallocate(phi_fmm, fx_fmm, fy_fmm)
  deallocate(phi_direct, fx_direct, fy_direct)

contains

  !============================================================================
  ! ПРЯМОЙ РАСЧЕТ φ, F
  !============================================================================
  subroutine compute_direct(N, xs, ys, qs, pot, phi, fx, fy)
    integer, intent(in) :: N
    real(dp), intent(in) :: xs(N), ys(N), qs(N)
    type(coulomb_2d_potential), intent(in) :: pot
    real(dp), intent(out) :: phi(N), fx(N), fy(N)

    integer :: i, j
    real(dp) :: dx, dy, r2, r, r_inv, r3_inv

    phi = 0.0_dp
    fx = 0.0_dp
    fy = 0.0_dp

    do i = 1, N
      do j = 1, N
        if (i == j) cycle

        dx = xs(i) - xs(j)
        dy = ys(i) - ys(j)
        r2 = dx*dx + dy*dy

        if (r2 < 1.0e-16_dp) cycle

        r = sqrt(r2)
        r_inv = 1.0_dp / r
        r3_inv = r_inv / r2

        phi(i) = phi(i) + qs(j) * r_inv
        fx(i) = fx(i) + qs(j) * dx * r3_inv
        fy(i) = fy(i) + qs(j) * dy * r3_inv
      end do
    end do

  end subroutine compute_direct

end program test_fmm_v2
