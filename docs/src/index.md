# Welcome to thermal-comforts's documentation!

```{toctree}
---
maxdepth: 1
```

The `thermal-comfort` package wraps a few common thermal-comfort functions from official
sources such as ISO-norms or VDI-Guidelines in python. The underlying functions are implemented in
fortran to achieve blazingly fast performance on large arrays. If possible, the original code was
reused and slightly modified to create a standardized interface for all function.

Tests are implemented for all functions using offical data often provided alongside a paper or the implementation itself.

## Installation

via https

```bash
pip install git+https://github.com/RUBclim/thermal-comfort
```

via ssh

```bash
pip install git+ssh://git@github.com/RUBclim/thermal-comfort
```

## Quick start

The termal-comfort package provides a limited set of commonly used functions.
Which work for scalar values, but are mainly optimized for large array calculation
of hundres of thousands of values.

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

The functions only accept 1-dimensional arrays, multi dimensional arrays must be reshaped before and after.

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

## API documentation

```{eval-rst}
.. automodule:: thermal_comfort
   :members:
   :undoc-members:
   :private-members:
```

## Indices and tables

- {ref}`genindex`
