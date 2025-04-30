[![ci](https://github.com/RUBclim/thermal-comfort/actions/workflows/ci.yml/badge.svg)](https://github.com/RUBclim/thermal-comfort/actions/workflows/ci.yml)
[![wheels](https://github.com/RUBclim/thermal-comfort/actions/workflows/wheels.yml/badge.svg)](https://github.com/RUBclim/thermal-comfort/actions/workflows/wheels.yml)
[![pre-commit](https://github.com/RUBclim/thermal-comfort/actions/workflows/pre-commit.yaml/badge.svg)](https://github.com/RUBclim/thermal-comfort/actions/workflows/pre-commit.yaml)

# thermal-comfort

## Installation

via pypi

```bash
pip install thermal-comfort
```

via https

```bash
pip install git+https://github.com/RUBclim/thermal-comfort
```

via ssh

```bash
pip install git+ssh://git@github.com/RUBclim/thermal-comfort
```

> [!NOTE]
> For this to work, you will have to have `git` and `gfortran` installed

For every release, pre-compiled ABI-3 wheels are provided under
[releases](https://github.com/RUBclim/thermal-comfort/releases)

## Run the tests

On version 3.12 for example...

```bash
tox -e py312
```

## Documentation

Docs can be found here: https://rubclim.github.io/thermal-comfort/.

## Quick start

The thermal-comfort package provides a limited set of commonly used functions. Which
work for scalar values, but are mainly optimized for large array calculation of hundreds
of thousands of values.

### scalars

```python
from thermal_comfort import utci_approx

utci_approx(ta=20.3, tmrt=50.9, v=2.7, rh=50.5)

```

### arrays

#### 1-dimensional arrays

```python
import numpy as np
from thermal_comfort import utci_approx

utci_approx(
    ta=np.array([20.3, 28.5]),
    tmrt=np.array([50.9, 70.3]),
    v=np.array([2.7, 1.9]),
    rh=np.array([50.5, 70.3]),
)

```

#### n-dimensional arrays

The functions only accept 1-dimensional arrays, multi dimensional arrays must be
reshaped before and after.

```python
import numpy as np
from thermal_comfort import utci_approx

# 2D arrays e.g. a raster
ta = np.array([[20.3, 28.5], [20.3, 28.5]])
tmrt = np.array([[50.9, 70.3], [50.9, 70.3]])
v = np.array([[2.7, 1.9], [2.7, 1.9]])
rh = np.array([[50.5, 70.3], [50.5, 70.3]])
# retrieve the initial shape
orig_shape = ta.shape

# reshape the array to be 1-dimensional
ta = np.ravel(ta)
tmrt = np.ravel(tmrt)
v = np.ravel(v)
rh = np.ravel(rh)

# calculate the UTCI along the 1-dimensional array
utci = utci_approx( ta=ta, tmrt=tmrt, v=v, rh=rh)

# restore the original shape
utci = utci.reshape(orig_shape)
```

### API

For a complete documentation look at the
[docs](https://rubclim.github.io/thermal-comfort/.)

#### Mean Radiant Temperature (MRT)

Calculate the mean radiant temperature based on DIN EN ISO 7726.

```python
mean_radiant_temp(ta, tg, v, d = 0.15, e = 0.95)
```

- `ta`: air temperature in °C
- `tg`: black globe temperature in °C

#### UTCI

Calculate the Universal Thermal Climate Index (UTCI)

```python
utci_approx(ta, tmrt, v, rh)
```

- `ta`: Air temperature in °C
- `tmrt`: Mean radiant temperature in °C
- `v`: Wind speed in m/s
- `rh`: Relative humidity in %

#### PET

Calculate the Physiological Equivalent Temperature (PET).

```python
pet_static(ta, tmrt, v, rh, p)
```

- `ta`: air temperature in °C
- `rh`: relative humidity in %
- `v`: wind speed in m/s
- `tmrt`: mean radiant temperature in °C
- `p`: atmospheric pressure in hPa

#### Heat Index

Calculate the heat index following Steadman R.G (1979) & Rothfusz L.P (1990).

```python
heat_index(ta, rh)
```

- `ta`: air temperature in °C
- `rh`: relative humidity in %

#### Extended Heat Index

Calculate the heat index following Steadman R.G (1979) & Rothfusz L.P (1990), but
extends the range following The National Weather Service Weather Predicion Center.

```python
heat_index_extended(ta, rh)
```

- `ta`: air temperature in °C
- `rh`: relative humidity in %

#### Wet Bulb Temperature (TWB)

Calculate the wet bulb temperature following the Stull (2011) equation

```python
wet_bulb_temp(ta, rh)
```

- `ta`: air temperature in °C
- `rh`: relative humidity in %

#### Saturation Vapor Pressure

##### Over water

```python
sat_vap_press_water(ta)
```

- `ta`: air temperature in °C

##### Over ice

```python
sat_vap_press_ice(ta)
```

- `ta`: air temperature in °C

#### Dew Point Temperature

```python
dew_point(ta, rh)
```

- `ta`: air temperature in °C
- `rh`: relative humidity in %

#### Absolute Humidity

```python
absolute_humidity(ta, rh)
```

- `ta`: air temperature in °C
- `rh`: relative humidity in %

#### Specific Humidity

```python
specific_humidity(ta, rh)
```

- `ta`: air temperature in °C
- `rh`: relative humidity in %

## Performance

The benchmark was ran using an array of 100,000 values and 4 threads
(`OMP_NUM_THREADS=4`). See the `benchmark` directory for more details on the benchmarks.

The hardware used is:

- 2x AMD EPYC 7702 64-Core Processor
- ubuntu 22.04

using an array of length 100,000 the following results were found:

### comparing to pythermalcomfort

| Benchmark         | pythermalcomfort ([bf9febd][1]) |     thermal-comfort      | thermal-comfort (unsafe) | thermal-comfort (Open MPI) | thermal-comfort (unsafe & Open MPI) |
| ----------------- | :-----------------------------: | :----------------------: | :----------------------: | :------------------------: | :---------------------------------: |
| tmrt scalar       |             21.0 us             |  20.3 us: 1.03x faster   |  2.21 us: 9.50x faster   |   23.3 us: 1.11x slower    |        4.38 us: 4.79x faster        |
| tmrt array        |             11.9 ms             |  5.77 ms: 2.06x faster   |  11.3 ms: 1.05x faster   |   3.86 ms: 3.07x faster    |        1.29 ms: 9.17x faster        |
| twb scalar        |             11.2 us             |  2.43 us: 4.59x faster   |  1.31 us: 8.52x faster   |   2.44 us: 4.57x faster    |        1.27 us: 8.80x faster        |
| twb array         |             8.86 ms             |  2.54 ms: 3.49x faster   |  2.52 ms: 3.52x faster   |   2.63 ms: 3.37x faster    |        1.73 ms: 5.13x faster        |
| heat index scalar |             34.9 us             |  2.38 us: 14.67x faster  |  1.24 us: 28.10x faster  |   5.06 us: 6.89x faster    |        3.68 us: 9.47x faster        |
| heat index array  |             4.06 ms             |  2.27 ms: 1.79x faster   |  2.44 ms: 1.66x faster   |   1.92 ms: 2.11x faster    |        859 us: 4.72x faster         |
| utci scalar       |             59.5 us             |  22.7 us: 2.62x faster   |  2.23 us: 26.66x faster  |   26.0 us: 2.29x faster    |       4.45 us: 13.36x faster        |
| utci array        |             37.1 ms             |  5.93 ms: 6.25x faster   |  11.4 ms: 3.25x faster   |   3.82 ms: 9.71x faster    |       1.24 ms: 29.84x faster        |
| pet scalar        |             6.04 ms             | 5.86 us: 1031.44x faster | 6.73 us: 897.57x faster  |  8.68 us: 696.14x faster   |       6.56 us: 922.14x faster       |
| pet array         |             189 sec             |  284 ms: 665.14x faster  |  804 ms: 234.82x faster  |  117 ms: 1611.19x faster   |       113 ms: 1669.13x faster       |
| Geometric mean    |              (ref)              |      10.00x faster       |      13.83x faster       |       10.45x faster        |            23.64x faster            |

### comparing to umep

| Benchmark      |   umep   |      umep (`njit`)      | pythermalcomfort ([bf9febd][1]) |     thermal-comfort      | thermal-comfort (unsafe) | thermal-comfort (Open MPI) | thermal-comfort (unsafe & Open MPI) |
| -------------- | :------: | :---------------------: | :-----------------------------: | :----------------------: | :----------------------: | :------------------------: | :---------------------------------: |
| utci scalar    | 44.5 us  |  6.56 us: 6.79x faster  |      59.5 us: 1.34x slower      |  22.7 us: 1.96x faster   |  2.23 us: 19.97x faster  |   26.0 us: 1.72x faster    |       4.45 us: 10.00x faster        |
| utci array     | 6.87 sec |  313 ms: 21.95x faster  |     37.1 ms: 185.23x faster     | 5.93 ms: 1157.28x faster | 11.4 ms: 602.16x faster  |  3.82 ms: 1798.79x faster  |      1.24 ms: 5526.68x faster       |
| pet scalar     |  388 us  | 3.17 us: 122.51x faster |     6.04 ms: 15.56x slower      |  5.86 us: 66.29x faster  |  6.73 us: 57.69x faster  |   8.68 us: 44.74x faster   |       6.56 us: 59.26x faster        |
| pet array      | 62.6 sec | 1.33 sec: 48.25x faster |      189 sec: 3.02x slower      |  284 ms: 220.43x faster  |  804 ms: 77.82x faster   |   117 ms: 533.95x faster   |       113 ms: 553.15x faster        |
| Geometric mean |  (ref)   |      35.61x faster      |          2.07x faster           |      53.17x faster       |      88.52x faster       |       51.68x faster        |           148.51x faster            |

While `njit` already gives a huge performance boost, the difference between umep
(`njit`) and thermal-comfort increases for larger arrays e.g. 1,000,000 values as shown
here:

| Benchmark      | umep (`njit`) |    thermal-comfort     |
| -------------- | :-----------: | :--------------------: |
| utci scalar    |    6.53 us    | 22.9 us: 3.51x slower  |
| utci array     |   3.09 sec    | 48.5 ms: 63.77x faster |
| pet scalar     |    3.15 us    | 5.84 us: 1.85x slower  |
| pet array      |   13.3 sec    | 2.80 sec: 4.75x faster |
| Geometric mean |     (ref)     |      2.61x faster      |

> [!CAUTION]
> If you're after the last bit of performance and don't care about input
> validation, you may use the underscored functions e.g. `_utci_approx` or `_pet_static`
> which fully avoid any computations in python. However, you will have to guarantee that
> all your arrays have the same length otherwise undefined behavior may happen. For
> performance reasons this package is not compiled using `-fcheck=bounds` compiler-flag.

[1]:
  https://github.com/CenterForTheBuiltEnvironment/pythermalcomfort/commit/bf9febdfb6244fff0fd9805c0ed1b41820504696
