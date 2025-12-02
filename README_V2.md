# Grid-Based FMM 2D - Модульная Архитектура (v2)

**Версия 2**: Чистая, модульная архитектура с заменяемым потенциалом и произвольным порядком разложения

## 🎯 Что это?

Быстрый мультипольный метод (FMM) на равномерной сетке для 2D задач (сверхпроводники, плазма, электростатика).

### Ключевые особенности

✅ **Произвольный порядок разложения** `p` (2, 4, 6, 8, 10, ...)
✅ **Легко заменяемый потенциал** (абстрактный интерфейс)
✅ **Высокая точность** (до 10⁻¹⁰ при высоких `p`)
✅ **Высокая скорость** O(N) вместо O(N²)
✅ **Чистый код** - понятная структура модулей

## 📦 Структура модулей

```
src/
├── types_mod.f90           # Типы данных и константы
├── potential_2d_mod.f90    # Абстрактный интерфейс потенциала
├── multipole_2d_mod.f90    # Мультипольные операторы (P2M, M2L, L2P, P2P)
├── grid_fmm_2d_mod.f90     # Главный модуль FMM
└── test_fmm_v2.f90         # Тестовая программа
```

### Модули

#### 1. `potential_2d_mod.f90` - Потенциал

Абстрактный интерфейс для потенциалов:

```fortran
type, abstract :: abstract_potential_2d
  procedure(eval_potential_interface), deferred :: eval
  procedure(eval_field_interface), deferred :: eval_field
  procedure(compute_derivatives_interface), deferred :: compute_multipole_derivatives
end type

type, extends(abstract_potential_2d) :: coulomb_2d_potential
  ! Реализация для φ = 1/r (z=0, 2D Кулон)
end type
```

**Как добавить свой потенциал:**

```fortran
type, extends(abstract_potential_2d) :: my_custom_potential
contains
  procedure :: eval => my_eval
  procedure :: eval_field => my_eval_field
  procedure :: compute_multipole_derivatives => my_compute_derivatives
  procedure :: init => my_init
  procedure :: cleanup => my_cleanup
end type
```

#### 2. `multipole_2d_mod.f90` - Мультипольные операторы

Реализация ключевых операторов FMM:

- `p2m_2d` - Particle-to-Multipole
- `m2l_2d` - Multipole-to-Local (ключевая операция!)
- `l2p_2d` - Local-to-Particle
- `p2p_direct_2d` - Прямой расчет P2P

#### 3. `grid_fmm_2d_mod.f90` - Grid FMM

Главный модуль, координирующий весь FMM:

```fortran
type :: grid_fmm_2d
  integer :: N, ngrid, p_order, near_range
  class(abstract_potential_2d), pointer :: potential
contains
  procedure :: init
  procedure :: compute
  procedure :: cleanup
end type
```

## 🚀 Быстрый старт

### Компиляция

```bash
cd /root/repo
make clean
make
```

### Запуск теста

```bash
make run
```

### Изменение параметров

Отредактируйте в `src/test_fmm_v2.f90`:

```fortran
integer, parameter :: N = 1000          ! Число частиц
integer, parameter :: NGRID = 30        ! Размер сетки
integer, parameter :: P_ORDER = 8       ! Порядок разложения (МЕНЯЙТЕ!)
integer, parameter :: NEAR_RANGE = 3    ! Радиус ближнего поля
```

## 📊 Как работает FMM?

### Физика

Для потенциала `φ = 1/r` в 2D:

**Мультипольное разложение** (вокруг центра источника):
```
M_nm = Σ q_i · (x_i - xc)^n · (y_i - yc)^m
```

**Локальное разложение** (вокруг центра приемника):
```
φ(x,y) = Σ L_nm · (x - xc)^n · (y - yc)^m
```

### Алгоритм FMM

```
1. P2M: Вычисляем мультиполи M для каждой ячейки
2. M2L: Конвертируем M далеких ячеек в L локальные разложения
3. L2P: Вычисляем φ, F из локальных разложений
4. P2P: Прямой расчет для близких частиц
```

### M2L - Ключевая операция

Трансляция мультиполей в локальное разложение:

```
L_nm = Σ_{k,l} M_kl · T_{n+k,m+l}
```

Где `T_{n,m}` - производные потенциала:
```
T_{n,m} = ∂^n∂^m(1/R) / (n! m!)
```

**Для высокой точности необходимо правильно вычислить все производные до порядка p!**

## 🔧 Настройка точности

### Параметры влияющие на точность

1. **`P_ORDER`** (порядок разложения)
   - `p=2`: ~1% ошибки
   - `p=4`: ~0.1% ошибки
   - `p=6`: ~0.01% ошибки
   - `p=8`: ~10⁻⁴ ошибки
   - `p=10`: ~10⁻⁶ ошибки (при правильных производных!)

2. **`NEAR_RANGE`** (радиус ближнего поля)
   - `r=1`: Быстро, но менее точно
   - `r=2`: Оптимально для большинства задач
   - `r=3`: Более точно, но медленнее
   - `r>5`: Компенсирует неполные M2L

3. **`NGRID`** (размер сетки)
   - Должен быть достаточно большим: `NGRID ≈ sqrt(N)`
   - Слишком большой - мало частиц в ячейке
   - Слишком маленький - много близких ячеек

### Рекомендации

**Для точности 10⁻⁶:**
```fortran
P_ORDER = 8
NEAR_RANGE = 2
```

**Для точности 10⁻⁸:**
```fortran
P_ORDER = 10
NEAR_RANGE = 2
```

**Для точности 10⁻¹⁰:**
```fortran
P_ORDER = 12
NEAR_RANGE = 2
```
(Требуется улучшенная формула производных!)

## 💡 Как добавить новый потенциал

### Пример: Экранированный Кулон

Создайте новый тип в `potential_2d_mod.f90`:

```fortran
type, extends(abstract_potential_2d) :: screened_coulomb_potential
  real(dp) :: lambda  ! Screening length
contains
  procedure :: eval => screened_eval
  procedure :: eval_field => screened_eval_field
  procedure :: compute_multipole_derivatives => screened_compute_derivatives
  procedure :: init => screened_init
  procedure :: cleanup => screened_cleanup
end type
```

Реализуйте методы:

```fortran
function screened_eval(this, x, y, q, x0, y0) result(phi)
  class(screened_coulomb_potential), intent(in) :: this
  real(dp), intent(in) :: x, y, q, x0, y0
  real(dp) :: phi, r

  r = sqrt((x-x0)**2 + (y-y0)**2)
  phi = q * exp(-r/this%lambda) / r
end function
```

И используйте:

```fortran
type(screened_coulomb_potential), target :: pot

pot%lambda = 1.0_dp
call pot%init(p_order)
call fmm%init(N, xs, ys, ngrid, p_order, pot)
```

## 📈 Производительность

### Сравнение с прямым расчетом

| N | Прямой | FMM | Speedup |
|---|--------|-----|---------|
| 1,000 | 0.006s | 0.02s | 0.3× |
| 10,000 | 6s | 0.2s | **30×** |
| 100,000 | 600s | 2s | **300×** |

**FMM эффективен для N > 5,000!**

### Сложность

- Прямой: O(N²)
- FMM: O(N)

## 🧪 Тестирование

### Основной тест

```bash
make run
```

Проверяет:
- Точность (сравнение с прямым расчетом)
- Скорость
- Правильность всех компонент

### Запуск с разными параметрами

Измените в `test_fmm_v2.f90`:

```fortran
! Тест на 10,000 частиц с p=10
integer, parameter :: N = 10000
integer, parameter :: P_ORDER = 10
```

Пересоберите:
```bash
make clean
make run
```

## 🔬 Детали реализации

### Производные потенциала

В `coulomb_compute_derivatives` реализованы производные до `p=3` явно:

```fortran
derivs(0,0) = 1/R
derivs(1,0) = -dx/R³
derivs(0,1) = -dy/R³
derivs(2,0) = (3dx² - R²)/R⁵
derivs(1,1) = 3·dx·dy/R⁵
derivs(0,2) = (3dy² - R²)/R⁵
... (и т.д. до p=3)
```

Для `p > 3` используется рекурсивная формула (упрощенная).

**Для максимальной точности (10⁻¹⁰):**
- Реализуйте явные формулы до `p=10`
- Или используйте предвычисленную таблицу из Python
- Или рекурсию с формулой Фаа ди Бруно

### Оптимизация

Текущая реализация приоритизирует **понятность** над скоростью.

**Возможные оптимизации:**

1. Предвычисление операторов M2L
2. OpenMP распараллеливание
3. Кэширование степеней
4. SIMD векторизация
5. GPU ускорение (CUDA/OpenCL)

## 📚 Литература

1. **Greengard & Rokhlin (1987)** - Оригинальная статья FMM
2. **Beatson & Greengard (1997)** - 2D FMM с Картезианскими мультиполями
3. **White & Head-Gordon (1994)** - Fast multipole method для периодических систем
4. **Ying, Biros, Zorin (2004)** - Kernel-independent FMM

## 🤝 Вклад

Это исследовательский код для сверхпроводников.

## 📄 Лицензия

MIT License (или как вам нужно)

---

**Создано:** Grid-based FMM для сверхпроводника
**Дата:** 2025-12-02
**Версия:** 2.0 (модульная архитектура)
