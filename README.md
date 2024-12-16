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

The vectorized fortran version is on average 4.2x faster than the `pythermalcomfort` version. An individual function call is 54.4x times faster.
