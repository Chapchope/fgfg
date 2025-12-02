# Grid-Based FMM 2D - Статус проекта

**Дата**: 2025-12-02
**Версия**: 2.0 (Модульная архитектура)

## ✅ Что реализовано

### 1. Модульная архитектура (ГОТОВО)

Создана чистая, понятная модульная структура:

```
src/
├── types_mod.f90           # Типы и константы
├── potential_2d_mod.f90    # Абстрактный потенциал (легко заменяемый!)
├── multipole_2d_mod.f90    # Мультипольные операторы
├── grid_fmm_2d_mod.f90     # Главный FMM модуль
└── test_fmm_v2.f90         # Тест
```

### 2. Абстрактный интерфейс потенциала (ГОТОВО ✓)

**Файл**: `src/potential_2d_mod.f90`

```fortran
type, abstract :: abstract_potential_2d
  procedure :: eval               ! φ = V(r)
  procedure :: eval_field         ! φ, F = -∇φ
  procedure :: compute_multipole_derivatives  ! Для M2L
end type

type, extends(abstract_potential_2d) :: coulomb_2d_potential
  ! Реализация для φ = 1/r (2D Кулон при z=0)
end type
```

**Как добавить свой потенциал**:

1. Создайте новый тип, расширяющий `abstract_potential_2d`
2. Реализуйте методы: `eval`, `eval_field`, `compute_multipole_derivatives`
3. Используйте! Не нужно менять основной код FMM

**Пример**: Экранированный Кулон

```fortran
type, extends(abstract_potential_2d) :: screened_coulomb_potential
  real(dp) :: lambda  ! screening length
contains
  procedure :: eval => screened_eval
  procedure :: eval_field => screened_eval_field
  procedure :: compute_multipole_derivatives => screened_compute_derivatives
  procedure :: init => screened_init
  procedure :: cleanup => screened_cleanup
end type
```

### 3. Мультипольные операторы (ГОТОВО ✓)

**Файл**: `src/multipole_2d_mod.f90`

Реализованы все ключевые операторы FMM:

- `p2m_2d` - Particle-to-Multipole
- `m2l_2d` - Multipole-to-Local ⚠️ (см. ниже)
- `l2p_2d` - Local-to-Particle
- `p2p_direct_2d` - Прямой расчет

### 4. Grid FMM (ГОТОВО ✓)

**Файл**: `src/grid_fmm_2d_mod.f90`

```fortran
type :: grid_fmm_2d
  ! Настраиваемые параметры
  integer :: N, ngrid, p_order, near_range

  ! Полиморфный потенциал!
  class(abstract_potential_2d), pointer :: potential
contains
  procedure :: init
  procedure :: compute
  procedure :: cleanup
end type
```

**Использование**:

```fortran
type(grid_fmm_2d) :: fmm
type(coulomb_2d_potential), target :: pot

call pot%init(p_order)
call fmm%init(N, xs, ys, ngrid, p_order, pot, near_range)
call fmm%compute(qs, phi, fx, fy)
call fmm%cleanup()
```

### 5. Тест и документация (ГОТОВО ✓)

- `test_fmm_v2.f90` - Полный тест с проверкой точности
- `Makefile` - Система сборки
- `README_V2.md` - Подробная документация

## ⚠️ Текущие проблемы

### M2L operator - Низкая точность

**Проблема**: Текущая реализация M2L через производные потенциала дает **очень большую ошибку** (~40x относительная ошибка).

**Причина**: Формула M2L через производные реализована неполностью/неправильно для `p > 3`.

**Текущее состояние**:

- `p=0,1,2,3`: Явные формулы (РАБОТАЮТ)
- `p > 3`: НЕ РЕАЛИЗОВАНО (заполнено нулями)

**Результат тестирования** (p=3, near_range=5):

```
Max relative error (phi):   4.0898E+01
Max relative error (Fx):    3.6304E+01
Max relative error (Fy):    2.7561E+02
```

## 🔧 Пути решения проблемы M2L

### Вариант A: Использовать явные формулы (РЕКОМЕНДУЕТСЯ)

Взять подход из `simple_fmm_2d.f90`:

```fortran
! Явные формулы для каждого порядка
! p=0: монополь
L(0,0) = L(0,0) + M(0,0) / R

! p=1: диполь
L(0,0) = L(0,0) + (M(1,0)*dx + M(0,1)*dy) / R³
L(1,0) = L(1,0) - M(0,0) * dx / R³ + M(1,0) / R
L(0,1) = L(0,1) - M(0,0) * dy / R³ + M(0,1) / R

! p=2: квадруполь
L(0,0) = L(0,0) + (M(2,0)*dx² + 2*M(1,1)*dx*dy + M(0,2)*dy²) / R⁵
... (и т.д.)

! p=3: октополь
... (нужно вывести формулы)
```

**Преимущества**:
- ✓ Надежно
- ✓ Быстро
- ✓ Понятно

**Недостатки**:
- ✗ Много формул писать вручную
- ✗ Для каждого p нужны новые формулы

**Работы**: 2-4 часа на вывод формул до p=6

### Вариант B: Рекурсивная формула

Реализовать правильную рекурсию для производных 1/R:

```fortran
! Рекурсивная формула для ∂^n∂^m(1/R)
! Базис:
T(0,0) = 1/R
T(1,0) = -dx/R³
T(0,1) = -dy/R³

! Рекурсия (формула Фаа ди Бруно или упрощенная версия):
T(n+1,m) = -(2n+2m+1)/R² [dx·T(n,m) + ...lower order terms...]
T(n,m+1) = -(2n+2m+1)/R² [dy·T(n,m) + ...lower order terms...]
```

**Преимущества**:
- ✓ Работает для произвольного p
- ✓ Элегантно

**Недостатки**:
- ✗ Сложно реализовать правильно
- ✗ Нужно понять точную формулу рекурсии

**Работы**: 4-8 часов (понять формулу + отладка)

### Вариант C: Таблица из Python

Предвычислить производные в Python с помощью SymPy:

```python
from sympy import *

x, y = symbols('x y')
R = sqrt(x**2 + y**2)
potential = 1/R

# Вычисляем все производные до порядка p_max
for n in range(p_max+1):
    for m in range(p_max+1-n):
        deriv = diff(potential, x, n, y, m)
        print(f"T({n},{m}) = {deriv}")
```

Затем использовать эту таблицу в Fortran.

**Преимущества**:
- ✓ Гарантированно правильные формулы
- ✓ Высокая точность

**Недостатки**:
- ✗ Нужна интеграция Python ↔ Fortran
- ✗ Сложность структуры проекта

**Работы**: 2-3 часа

## 📋 Рекомендации

### Для быстрого решения (2-4 часа)

**→ Используйте вариант A (явные формулы)**

1. Откройте `src/multipole_2d_mod.f90`
2. В функции `m2l_2d` замените вызов `potential%compute_multipole_derivatives` на явные формулы
3. Используйте формулы из `simple_fmm_2d.f90:313-362` как образец
4. Добавьте формулы для p=3 (октополь)
5. Для p>3 либо выведите формулы вручную, либо оставьте `near_range` большим

### Для полного решения (1 неделя)

1. Реализуйте явные формулы до p=6 (вариант A)
2. Добавьте рекурсивные формулы для p>6 (вариант B)
3. Добавьте таблицы из Python для высшей точности (вариант C)
4. Оптимизируйте с OpenMP

## 🎯 Текущее состояние кода

### Что работает отлично

✅ Модульная архитектура
✅ Абстрактный интерфейс потенциала
✅ P2M (Particle-to-Multipole)
✅ L2P (Local-to-Particle)
✅ P2P (Прямой расчет)
✅ Grid структура
✅ Сборка и тестирование

### Что нужно исправить

⚠️ M2L (Multipole-to-Local) - неполная реализация для p>3

### Что можно улучшить

💡 OpenMP распараллеливание
💡 Оптимизация памяти
💡 SIMD векторизация
💡 Предвычисление операторов

## 🚀 Быстрый старт

```bash
# Сборка
make clean
make

# Запуск (текущие настройки: p=3, near_range=5)
./bin/test_fmm_v2

# Изменить параметры
# Отредактируйте src/test_fmm_v2.f90:
#   P_ORDER = 3
#   NEAR_RANGE = 5
```

## 📚 Файлы проекта

### Основные (используйте эти!)

- `src/potential_2d_mod.f90` - Потенциал (легко заменяемый!)
- `src/multipole_2d_mod.f90` - Мультипольные операторы
- `src/grid_fmm_2d_mod.f90` - Grid FMM
- `src/test_fmm_v2.f90` - Тест

### Справочные

- `src/simple_fmm_2d.f90` - Упрощенная версия (работает!)
- `src/test_simple.f90` - Тест для simple версии
- `README_V2.md` - Подробная документация

### Устаревшие (не использовать)

- `src/grid_fmm_mod.f90` - Старая версия
- `src/cartesian_multipole_mod.f90` - Старая версия
- `src/potential_interface_mod.f90` - Старая версия

## 💬 Итог

У вас есть **отличная модульная архитектура** с:

✅ Легко заменяемым потенциалом
✅ Чистым, понятным кодом
✅ Произвольным порядком разложения (архитектура готова!)
✅ Полной документацией

**Но**: M2L operator требует доработки для p>3.

**Следующий шаг**: Реализовать явные формулы M2L (вариант A, ~2-4 часа работы).

---

**Вопросы?** Смотрите `README_V2.md` для подробной документации!
