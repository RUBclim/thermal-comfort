[![ci](https://github.com/RUBclim/thermal-comfort/actions/workflows/ci.yml/badge.svg)](https://github.com/RUBclim/thermal-comfort/actions/workflows/ci.yml)
[![pre-commit](https://github.com/RUBclim/thermal-comfort/actions/workflows/pre-commit.yaml/badge.svg)](https://github.com/RUBclim/thermal-comfort/actions/workflows/pre-commit.yaml)

# thermal-comfort

## Installation

via https

```bash
pip install git+https://github.com/RUBclim/thermal-comfort
```

via ssh

```bash
pip install git+ssh://git@github.com/RUBclim/thermal-comfort
```

## Run the tests

On version 3.12 for example...

```bash
tox -e py312
```

## performance

The array benchmarks are run on an array 1D-array of length 5000. With longer arrays the
difference improves. Due to the PET-calculation being so slow in `pythermalcomfort` we
cannot use larger arrays.

| Benchmark      | pythermalcomfort |     thermalcomfort      |
| -------------- | :--------------: | :---------------------: |
| tmrt scalar    |     18.9 us      |  3.64 us: 5.21x faster  |
| tmrt array     |      513 us      |  987 us: 1.92x slower   |
| twb scalar     |     10.5 us      |  2.03 us: 5.19x faster  |
| twb array      |      397 us      |  192 us: 2.07x faster   |
| utci scalar    |     31.7 us      |  3.58 us: 8.87x faster  |
| utci array     |     2.17 ms      |  746 us: 2.90x faster   |
| pet scalar     |     4.89 ms      | 8.76 us: 558.34x faster |
| pet array      |     16.4 sec     | 84.3 ms: 195.00x faster |
| Geometric mean |      (ref)       |      9.75x faster       |

using an array of length 100,000

| Benchmark      | pythermalcomfort_large_array | thermalcomfort_large_array |
| -------------- | :--------------------------: | :------------------------: |
| tmrt scalar    |           18.8 us            |   2.34 us: 8.05x faster    |
| tmrt array     |           11.8 ms            |   15.2 ms: 1.28x slower    |
| twb scalar     |           11.9 us            |   1.26 us: 9.40x faster    |
| twb array      |           9.81 ms            |   3.22 ms: 3.04x faster    |
| utci scalar    |           33.0 us            |   2.19 us: 15.07x faster   |
| utci array     |           46.6 ms            |   10.9 ms: 4.26x faster    |
| pet scalar     |           4.72 ms            |  5.57 us: 848.07x faster   |
| pet array      |                              |    1.02 sec +- 0.02 sec    |
| Geometric mean |            (ref)             |        9.97x faster        |
