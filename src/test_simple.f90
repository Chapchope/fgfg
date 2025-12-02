program test_simple
  use types_mod
  use simple_fmm_2d
  implicit none

  integer, parameter :: N = 1000
  real(dp) :: xs(N), ys(N), qs(N)
  real(dp) :: phi_fmm(N), fx_fmm(N), fy_fmm(N)
  real(dp) :: phi_dir(N), fx_dir(N), fy_dir(N)
  real(dp) :: t0, t1, t_fmm, t_dir
  real(dp) :: err_phi, err_f, max_phi, max_f
  integer :: i

  print *, "========================================"
  print *, "  Test Simple FMM 2D"
  print *, "========================================"
  print *, ""

  ! Инициализация случайных частиц
  call random_seed()
  call random_number(xs)
  call random_number(ys)

  xs = (xs - 0.5_dp) * 2.0_dp
  ys = (ys - 0.5_dp) * 2.0_dp
  qs = 1.0_dp

  print '(A,I0)', "Particles: ", N
  print *, ""

  ! FMM расчет
  print *, "Running FMM..."
  call cpu_time(t0)
  call fmm_2d_simple(N, xs, ys, qs, phi_fmm, fx_fmm, fy_fmm)
  call cpu_time(t1)
  t_fmm = t1 - t0
  print '(A,F8.4,A)', "FMM time: ", t_fmm, " s"
  print *, ""

  ! Прямой расчет (для проверки)
  print *, "Running direct sum (for validation)..."
  call cpu_time(t0)
  call direct_sum(N, xs, ys, qs, phi_dir, fx_dir, fy_dir)
  call cpu_time(t1)
  t_dir = t1 - t0
  print '(A,F8.4,A)', "Direct time: ", t_dir, " s"
  print '(A,F6.1,A)', "Speedup: ", t_dir/t_fmm, "x"
  print *, ""

  ! Сравнение точности
  print *, "========================================"
  print *, "ACCURACY CHECK"
  print *, "========================================"
  print *, ""

  ! Найдем максимальные ошибки
  max_phi = 0.0_dp
  max_f = 0.0_dp

  do i = 1, N
    err_phi = abs(phi_fmm(i) - phi_dir(i)) / max(abs(phi_dir(i)), 1.0_dp)
    err_f = sqrt((fx_fmm(i)-fx_dir(i))**2 + (fy_fmm(i)-fy_dir(i))**2) / &
            max(sqrt(fx_dir(i)**2 + fy_dir(i)**2), 1.0_dp)

    max_phi = max(max_phi, err_phi)
    max_f = max(max_f, err_f)
  end do

  print '(A,ES12.4)', "Max relative error (phi): ", max_phi
  print '(A,ES12.4)', "Max relative error (F):   ", max_f
  print *, ""

  ! Показываем примеры
  print *, "Sample comparison (first 5 particles):"
  print *, "   i  |   phi_FMM   |  phi_Direct | Rel. Error"
  print *, "------|-------------|-------------|------------"
  do i = 1, min(5, N)
    err_phi = abs(phi_fmm(i) - phi_dir(i)) / max(abs(phi_dir(i)), 1.0_dp)
    print '(I6," | ",F11.4," | ",F11.4," | ",ES10.2)', &
          i, phi_fmm(i), phi_dir(i), err_phi
  end do
  print *, ""

  print *, "========================================"
  if (max_phi < 1e-4_dp) then
    print *, "✓ Test PASSED! Accuracy is good."
  else if (max_phi < 1e-2_dp) then
    print *, "⚠ Test OK, but accuracy could be better."
  else
    print *, "✗ Test FAILED! Accuracy is poor."
  end if
  print *, "========================================"

contains

  subroutine direct_sum(N, xs, ys, qs, phi, fx, fy)
    integer, intent(in) :: N
    real(dp), intent(in) :: xs(N), ys(N), qs(N)
    real(dp), intent(out) :: phi(N), fx(N), fy(N)

    integer :: i, j
    real(dp) :: rx, ry, r2, r, r_inv, r3_inv

    phi = 0.0_dp
    fx = 0.0_dp
    fy = 0.0_dp

    do i = 1, N
      do j = 1, N
        if (i == j) cycle

        rx = xs(i) - xs(j)
        ry = ys(i) - ys(j)
        r2 = rx*rx + ry*ry

        if (r2 < 1e-16_dp) cycle

        r = sqrt(r2)
        r_inv = 1.0_dp / r
        r3_inv = r_inv / r2

        phi(i) = phi(i) + qs(j) * r_inv
        fx(i) = fx(i) + qs(j) * rx * r3_inv
        fy(i) = fy(i) + qs(j) * ry * r3_inv
      end do
    end do

  end subroutine direct_sum

end program test_simple
