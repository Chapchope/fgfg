program test_accurate
  !============================================================================
  ! ТЕСТ: Точный Grid-based FMM 2D
  !
  ! Проверяем точность и скорость
  !============================================================================
  use types_mod
  use grid_fmm_accurate
  implicit none

  integer, parameter :: N = 1000       ! Число частиц
  integer, parameter :: P_ORDER = 2    ! Порядок разложения (M2L полный до p=2!)
  integer, parameter :: NGRID = 30     ! Размер сетки
  integer, parameter :: NEAR_RANGE = 9 ! Радиус ближнего поля (большой!)

  real(dp), allocatable :: xs(:), ys(:), qs(:)
  real(dp), allocatable :: phi_fmm(:), fx_fmm(:), fy_fmm(:)
  real(dp), allocatable :: phi_direct(:), fx_direct(:), fy_direct(:)

  real(dp) :: t1, t2, time_fmm, time_direct
  real(dp) :: err_phi, err_fx, err_fy
  integer :: i

  print '(A)', ""
  print '(A)', "============================================"
  print '(A)', "  TEST: Accurate Grid-based FMM 2D"
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

  xs = xs * 10.0_dp - 5.0_dp
  ys = ys * 10.0_dp - 5.0_dp
  qs = (qs - 0.5_dp) * 2.0_dp

  print '(A,I0,A)', "Generated ", N, " random particles"
  print '(A)', ""

  !============================================================================
  ! 1. ПРЯМОЙ РАСЧЕТ
  !============================================================================
  print '(A)', "Computing direct (exact) solution..."

  call cpu_time(t1)
  call compute_direct(N, xs, ys, qs, phi_direct, fx_direct, fy_direct)
  call cpu_time(t2)
  time_direct = t2 - t1

  print '(A,F10.6,A)', "Direct time: ", time_direct, " s"
  print '(A)', ""

  !============================================================================
  ! 2. FMM РАСЧЕТ
  !============================================================================
  print '(A)', "Computing FMM solution..."

  call cpu_time(t1)
  call fmm_compute(N, xs, ys, qs, P_ORDER, NGRID, NEAR_RANGE, phi_fmm, fx_fmm, fy_fmm)
  call cpu_time(t2)
  time_fmm = t2 - t1

  print '(A,F10.6,A)', "FMM time:    ", time_fmm, " s"
  print '(A)', ""

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
    print '(A)', "✓ EXCELLENT! Error < 10^-6"
  else if (err_phi < 1.0e-4_dp) then
    print '(A)', "✓ GOOD! Error < 10^-4"
  else if (err_phi < 1.0e-2_dp) then
    print '(A)', "✓ ACCEPTABLE! Error < 1%"
  else
    print '(A)', "○ Moderate accuracy. Try increasing P_ORDER or NEAR_RANGE"
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
    print '(A)', "○ For N=1000, direct may be competitive."
    print '(A)', "  FMM advantage appears at N=10,000+"
  end if
  print '(A)', ""

  print '(A)', "============================================"
  print '(A)', ""

  deallocate(xs, ys, qs)
  deallocate(phi_fmm, fx_fmm, fy_fmm)
  deallocate(phi_direct, fx_direct, fy_direct)

contains

  !============================================================================
  ! ПРЯМОЙ РАСЧЕТ
  !============================================================================
  subroutine compute_direct(N, xs, ys, qs, phi, fx, fy)
    integer, intent(in) :: N
    real(dp), intent(in) :: xs(N), ys(N), qs(N)
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

end program test_accurate
