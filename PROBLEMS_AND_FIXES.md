# Проблемы с точностью и как их исправить

## Проблема №1: Неправильная реализация M2L

**Что не так в текущей версии:**

В файле `src/cartesian_multipole_mod.f90` функция `m2l_translate_2d` использует неправильную формулу для трансляции мультипольных моментов в локальное разложение.

**Корень проблемы:**

Для потенциала Кулона 1/r в 2D, M2L трансляция НЕ может быть сделана через простую комплексную конволюцию. Нужна специальная формула, учитывающая:

1. Геометрию 2D (не 3D!)
2. Правильное разложение Тейлора для 1/r
3. Операторы трансляции для каждого порядка отдельно

**Правильная формула** (упрощенная):

```
Для R = sqrt(dx² + dy²) (расстояние между центрами ячеек)

L₀₀ += M₀₀ / R
L₀₀ += (M₁₀·dx + M₀₁·dy) / R³
L₁₀ += -M₀₀·dx/R³ + M₁₀/R
L₀₁ += -M₀₀·dy/R³ + M₀₁/R
...
```

## Проблема №2: Сложная структура кода

**Что не так:**

- Слишком много уровней абстракции (potential_interface, abstract types)
- Неочевидна физика происходящего
- Трудно отладить

**Решение:**

Создать простой модуль `simple_fmm_2d.f90` где:
- Все понятно и прозрачно
- Каждая функция делает одно и хорошо документирована
- Физика очевидна из кода

## Как исправить

### Вариант 1: Исправить существующий код

1. Исправить `m2l_translate_2d` в `cartesian_multipole_mod.f90`:

```fortran
subroutine m2l_translate_2d_correct(M, dx, dy, R, L)
  real(dp), intent(in) :: M(0:p_order, 0:p_order)
  real(dp), intent(in) :: dx, dy, R
  real(dp), intent(inout) :: L(0:p_order, 0:p_order)

  real(dp) :: R_inv, R2, R3, R5, R7

  R2 = R * R
  R_inv = 1.0_dp / R
  R3 = R_inv / R2
  R5 = R3 / R2
  R7 = R5 / R2

  ! Monopole contribution
  L(0,0) = L(0,0) + M(0,0) * R_inv

  ! Dipole contribution to monopole
  L(0,0) = L(0,0) + (M(1,0)*dx + M(0,1)*dy) * R3

  ! Dipole contribution to dipole
  L(1,0) = L(1,0) - M(0,0) * dx * R3 + M(1,0) * R_inv
  L(0,1) = L(0,1) - M(0,0) * dy * R3 + M(0,1) * R_inv

  ! Quadrupole contributions
  if (p_order >= 2) then
    ! To monopole
    L(0,0) = L(0,0) + (M(2,0)*dx*dx + 2.0_dp*M(1,1)*dx*dy + M(0,2)*dy*dy) * R5

    ! To dipole
    L(1,0) = L(1,0) + (M(1,0)*dx + M(0,1)*dy) * (-dx) * 3.0_dp * R5 + M(2,0) * R3
    L(0,1) = L(0,1) + (M(1,0)*dx + M(0,1)*dy) * (-dy) * 3.0_dp * R5 + M(0,2) * R3

    ! To quadrupole
    L(2,0) = L(2,0) + M(0,0) * (3.0_dp*dx*dx - R2) * R5 + M(2,0) * R_inv
    L(1,1) = L(1,1) + M(0,0) * (3.0_dp*dx*dy) * R5 + M(1,1) * R_inv
    L(0,2) = L(0,2) + M(0,0) * (3.0_dp*dy*dy - R2) * R5 + M(0,2) * R_inv
  end if

end subroutine
```

2. Увеличить `near_range` до 2 или 3

3. Увеличить `p_order` до 8-10

### Вариант 2: Использовать простую версию

Я начал создавать `simple_fmm_2d.f90` с правильной реализацией, но там была ошибка компиляции.

**Чтобы исправить:**

1. Заменить все переменные цикла `n`, `m` на `nn`, `mm`
2. Переименовать параметр `cell` в `grid_cell` чтобы избежать конфликта с типом

## Ожидаемая точность

После исправления:

| p_order | near_range | Ожидаемая точность |
|---------|------------|-------------------|
| 4 | 2 | 10⁻⁴ |
| 6 | 2 | 10⁻⁶ |
| 8 | 2 | 10⁻⁸ |
| 10 | 3 | 10⁻¹⁰ |

## Тестирование

После исправления запустите:

```bash
make clean
make
./bin/test_grid_fmm
```

Проверьте что:
- Max relative error < 10⁻⁶ для p=6
- RMS error < 10⁻⁷ для p=6
- Speedup > 10× для N=10000

## Дополнительные улучшения

1. **Adaptive near-range**: Автоматически подбирать near_range в зависимости от p_order
2. **Precompute operators**: Предвычислить M2L операторы для всех смещений
3. **Vectorization**: Использовать SIMD инструкции для P2P
4. **OpenMP**: Распараллелить M2L и P2P циклы

## Литература

Для правильной реализации M2L в 2D см:

1. Beatson & Greengard (1997) - "A short course on fast multipole methods"
   Раздел про 2D Laplace equation

2. Carrier, Greengard & Rokhlin (1988) - "A fast adaptive multipole algorithm"
   Формулы трансляции для картезианских координат

3. White & Head-Gordon (1994) - "Derivation and efficient implementation of the fast multipole method"
   Явные формулы для низких порядков
