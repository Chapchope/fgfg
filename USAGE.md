# Руководство по использованию Grid-Based FMM

## Содержание

1. [Быстрый старт](#быстрый-старт)
2. [Настройка параметров](#настройка-параметров)
3. [Примеры использования](#примеры-использования)
4. [Создание нового потенциала](#создание-нового-потенциала)
5. [Оптимизация производительности](#оптимизация-производительности)
6. [Решение проблем](#решение-проблем)

## Быстрый старт

### 1. Сборка

```bash
make clean
make
```

### 2. Запуск демо

```bash
# Быстрая демонстрация (N=1000)
./bin/demo_quick

# Полный тест с валидацией (N=10000, медленно)
# ./bin/test_grid_fmm
```

### 3. Визуализация (опционально)

```bash
python3 visualize.py --file frame_00000.csv --output particles.png
```

## Настройка параметров

### Основные параметры FMM

```fortran
type(grid_fmm_2d_params) :: params

params%p_order = 6        ! Порядок мультипольного разложения
params%ngrid_x = 20       ! Количество ячеек по X
params%ngrid_y = 20       ! Количество ячеек по Y
params%near_range = 2     ! Радиус ближнего поля (в ячейках)
params%use_precompute = .true.  ! Предвычисление факториалов
params%box_margin = 0.1_dp      ! Отступ от границ
```

### Выбор порядка разложения

Порядок `p_order` определяет точность:

- **p = 4**: Грубая точность (~10^-4), быстро
- **p = 6**: Стандартная точность (~10^-6), рекомендуется для начала
- **p = 8**: Высокая точность (~10^-8)
- **p = 10**: Очень высокая точность (~10^-10)
- **p = 12+**: Экстремальная точность, медленно

**Рекомендация**: Начните с `p = 6`, затем увеличивайте если нужна большая точность.

### Выбор размера сетки

Эмпирическое правило: `ngrid ≈ sqrt(N)`

Примеры:
```
N = 1,000   → ngrid = 30
N = 10,000  → ngrid = 100
N = 100,000 → ngrid = 300
```

Можно использовать неквадратную сетку если частицы неравномерно распределены:
```fortran
params%ngrid_x = 50
params%ngrid_y = 100  ! Если частицы вытянуты по Y
```

### Выбор радиуса ближнего поля

- `near_range = 1`: Минимум прямых взаимодействий, требует высокий `p_order`
- `near_range = 2`: **Оптимальный баланс (рекомендуется)**
- `near_range = 3`: Высокая точность, больше прямых взаимодействий

## Примеры использования

### Пример 1: Минимальная программа

```fortran
program minimal_fmm
  use types_mod
  use potential_interface_mod
  use grid_fmm_mod
  implicit none

  integer, parameter :: N = 1000
  real(dp) :: xs(N), ys(N), qs(N)
  real(dp) :: phi(N), fx(N), fy(N)
  type(coulomb_2d_potential) :: pot
  type(grid_fmm_2d_params) :: par

  ! Инициализация частиц
  call random_number(xs)
  call random_number(ys)
  xs = xs - 0.5_dp
  ys = ys - 0.5_dp
  qs = 1.0_dp

  ! Настройка
  par%p_order = 6
  par%ngrid_x = 30
  par%ngrid_y = 30
  par%near_range = 2
  par%use_precompute = .true.

  ! Создание потенциала
  pot = create_coulomb_2d_potential()

  ! Вычисление
  call grid_fmm_2d_compute(N, xs, ys, qs, pot, par, phi, fx, fy)

  ! Использовать phi, fx, fy...
  print *, "Потенциал в точке 1:", phi(1)
  print *, "Сила в точке 1:", fx(1), fy(1)

end program minimal_fmm
```

### Пример 2: MD симуляция

```fortran
program md_simulation
  use types_mod
  use potential_interface_mod
  use grid_fmm_mod
  use md_utils_mod
  implicit none

  integer, parameter :: N = 5000
  integer, parameter :: nsteps = 100
  real(dp), parameter :: dt = 0.001_dp

  real(dp) :: xs(N), ys(N), qs(N)
  real(dp) :: vxs(N), vys(N)
  real(dp) :: phi(N), fx(N), fy(N)
  type(coulomb_2d_potential) :: pot
  type(grid_fmm_2d_params) :: par
  integer :: step

  ! Инициализация
  call random_number(xs)
  call random_number(ys)
  xs = xs - 0.5_dp
  ys = ys - 0.5_dp
  qs = 1.0_dp
  vxs = 0.0_dp
  vys = 0.0_dp

  ! Настройка FMM
  par%p_order = 6
  par%ngrid_x = 70
  par%ngrid_y = 70
  par%near_range = 2
  par%use_precompute = .true.

  pot = create_coulomb_2d_potential()

  ! MD цикл
  do step = 1, nsteps
    ! Вычислить силы
    call grid_fmm_2d_compute(N, xs, ys, qs, pot, par, phi, fx, fy)

    ! Интегрировать
    call md_step_euler(N, xs, ys, vxs, vys, fx, fy, dt)

    ! Граничные условия
    call apply_reflect_bc(N, xs, ys, vxs, vys, 1.0_dp)

    ! Сохранить кадр
    if (mod(step, 10) == 0) then
      call write_frame_csv(step, N, xs, ys)
    end if
  end do

end program md_simulation
```

### Пример 3: Исследование точности

```fortran
program accuracy_study
  use types_mod
  use potential_interface_mod
  use grid_fmm_mod
  implicit none

  integer, parameter :: N = 1000
  integer :: p
  real(dp) :: xs(N), ys(N), qs(N)
  real(dp) :: phi_fmm(N), fx_fmm(N), fy_fmm(N)
  real(dp) :: phi_ref(N), fx_ref(N), fy_ref(N)
  real(dp) :: error
  type(coulomb_2d_potential) :: pot
  type(grid_fmm_2d_params) :: par

  ! Инициализация
  call random_number(xs)
  call random_number(ys)
  xs = xs - 0.5_dp
  ys = ys - 0.5_dp
  qs = 1.0_dp

  pot = create_coulomb_2d_potential()

  ! Эталонное решение (прямой расчет)
  call direct_sum(N, xs, ys, qs, phi_ref, fx_ref, fy_ref)

  ! Тестируем разные порядки
  print *, "p_order | Max Error (phi) | RMS Error (phi)"
  print *, "--------|-----------------|----------------"

  do p = 2, 12, 2
    par%p_order = p
    par%ngrid_x = 30
    par%ngrid_y = 30
    par%near_range = 2
    par%use_precompute = .true.

    call grid_fmm_2d_compute(N, xs, ys, qs, pot, par, phi_fmm, fx_fmm, fy_fmm)

    error = maxval(abs(phi_fmm - phi_ref))
    print '(I7," | ",ES14.6," | ",ES14.6)', p, error, &
          sqrt(sum((phi_fmm - phi_ref)**2) / real(N, dp))
  end do

contains

  subroutine direct_sum(N, xs, ys, qs, phi, fx, fy)
    ! ... реализация прямого расчета ...
  end subroutine

end program accuracy_study
```

## Создание нового потенциала

### Шаблон для нового потенциала

```fortran
module yukawa_potential_mod
  use types_mod
  use potential_interface_mod
  implicit none

  ! Yukawa потенциал: phi = q * exp(-kappa*r) / r
  type, extends(potential_type) :: yukawa_potential
     real(dp) :: kappa = 1.0_dp  ! Screening parameter
   contains
     procedure :: value => yukawa_value
     procedure :: force => yukawa_force
     procedure :: moment_coefficient => yukawa_moment
  end type yukawa_potential

contains

  pure function yukawa_value(this, r, q) result(phi)
    class(yukawa_potential), intent(in) :: this
    real(dp), intent(in) :: r, q
    real(dp) :: phi

    if (r < TINY) then
      phi = 0.0_dp
    else
      phi = q * exp(-this%kappa * r) / r
    end if
  end function yukawa_value

  pure subroutine yukawa_force(this, rx, ry, q, fx, fy)
    class(yukawa_potential), intent(in) :: this
    real(dp), intent(in) :: rx, ry, q
    real(dp), intent(out) :: fx, fy

    real(dp) :: r, r2, r3, factor

    r2 = rx*rx + ry*ry

    if (r2 < TINY*TINY) then
      fx = 0.0_dp
      fy = 0.0_dp
      return
    end if

    r = sqrt(r2)
    r3 = r * r2

    ! F = -grad(phi) = q * exp(-k*r) * (1 + k*r) * r / r^3
    factor = q * exp(-this%kappa * r) * (1.0_dp + this%kappa * r) / r3

    fx = factor * rx
    fy = factor * ry
  end subroutine yukawa_force

  pure function yukawa_moment(this, n, m, r) result(coeff)
    class(yukawa_potential), intent(in) :: this
    integer, intent(in) :: n, m
    real(dp), intent(in) :: r
    real(dp) :: coeff

    ! Для Yukawa: требуется модифицированное разложение
    ! Здесь упрощенная версия
    if (r < TINY) then
      if (n == 0 .and. m == 0) then
        coeff = 1.0_dp
      else
        coeff = 0.0_dp
      end if
    else
      coeff = r**n * exp(-this%kappa * r)
    end if
  end function yukawa_moment

end module yukawa_potential_mod
```

### Использование нового потенциала

```fortran
program test_yukawa
  use yukawa_potential_mod
  use grid_fmm_mod
  implicit none

  type(yukawa_potential) :: pot
  ! ... остальная настройка ...

  pot%kappa = 2.0_dp  ! Установить параметр screening

  call grid_fmm_2d_compute(N, xs, ys, qs, pot, par, phi, fx, fy)
end program test_yukawa
```

## Оптимизация производительности

### 1. Выбор компилятора

```bash
# GNU Fortran (хорошо для разработки)
make FC=gfortran

# Intel Fortran (обычно быстрее на 20-40%)
make FC=ifort
```

### 2. Оптимальные параметры для разных N

| N | p_order | ngrid | near_range | Время (ориентировочно) |
|---|---------|-------|------------|------------------------|
| 1,000 | 6 | 30 | 2 | 0.01 с |
| 10,000 | 6 | 100 | 2 | 0.1 с |
| 100,000 | 6 | 300 | 2 | 1-2 с |
| 1,000,000 | 6 | 1000 | 2 | 20-30 с |

### 3. Профилирование

Используйте `gprof` для поиска узких мест:

```bash
# Компиляция с профилированием
make clean
make FFLAGS="-pg -O2"

# Запуск
./bin/test_grid_fmm

# Анализ
gprof bin/test_grid_fmm gmon.out > analysis.txt
```

### 4. Советы по оптимизации

1. **Предвычисление**: Всегда используйте `use_precompute = .true.`

2. **Размер сетки**: Не делайте слишком большую сетку
   - Оптимум: `ngrid ≈ sqrt(N)`
   - Слишком большая → много M2L операций
   - Слишком маленькая → много P2P операций

3. **Порядок разложения**: Не используйте выше чем нужно
   - p = 6 обычно достаточно для большинства приложений
   - p > 10 редко необходимо

4. **Near-field range**: near_range = 2 обычно оптимален

## Решение проблем

### Проблема: Низкая точность

**Решение:**
1. Увеличьте `p_order` (попробуйте 8 или 10)
2. Увеличьте `near_range` до 3
3. Проверьте что `use_precompute = .true.`

### Проблема: Медленное выполнение

**Решение:**
1. Уменьшите `p_order` (попробуйте 6 или 4)
2. Оптимизируйте `ngrid` (должен быть примерно sqrt(N))
3. Используйте Intel Fortran вместо GNU
4. Проверьте что компилятор использует флаги оптимизации

### Проблема: Численная нестабильность

**Решение:**
1. Проверьте что частицы не слишком близко друг к другу
2. Используйте double precision (по умолчанию в этой реализации)
3. Уменьшите `p_order` если он слишком высокий (>12)

### Проблема: Сегфолты или ошибки памяти

**Решение:**
1. Проверьте размеры массивов
2. Скомпилируйте в debug режиме: `make debug`
3. Используйте valgrind: `valgrind ./bin/test_grid_fmm`

## Дополнительные ресурсы

- README.md - Общая документация
- Исходный код в `src/` - Хорошо документирован
- Примеры программ: `src/test_grid_fmm.f90`, `src/demo_quick.f90`

## Обратная связь

Если у вас есть вопросы или предложения, пожалуйста создайте issue в репозитории.
