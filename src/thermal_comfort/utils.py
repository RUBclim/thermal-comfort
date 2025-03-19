from collections.abc import Iterable
from typing import overload
from typing import TypeVar
from typing import Union

import numpy as np
import numpy.typing as npt

from ._thermal_comfort import thermal_comfort_mod


T = TypeVar('T', bound=Union[np.floating, np.integer])


# autopep8: off
@overload
def mean_radiant_temp(
        ta: float,
        tg: float,
        v: float,
        d: float = 0.15,
        e: float = 0.95,
) -> float: ...


@overload
def mean_radiant_temp(
        ta: Iterable[float],
        tg: Iterable[float],
        v: Iterable[float],
        d: Union[Iterable[float], float] = 0.15,
        e: Union[Iterable[float], float] = 0.95,
) -> npt.NDArray[T]: ...


@overload
def mean_radiant_temp(
        ta: npt.NDArray[T],
        tg: npt.NDArray[T],
        v: npt.NDArray[T],
        d: Union[npt.NDArray[T], float] = 0.15,
        e: Union[npt.NDArray[T], float] = 0.95,
) -> npt.NDArray[T]: ...
# autopep8: on


def mean_radiant_temp(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
        tg: Union[npt.NDArray[T], float, Iterable[float]],
        v: Union[npt.NDArray[T], float, Iterable[float]],
        d: Union[npt.NDArray[T], float, Iterable[float]] = 0.15,
        e: Union[npt.NDArray[T], float, Iterable[float]] = 0.95,
) -> Union[npt.NDArray[T], float]:
    """
    Calculate the mean radiant temperature based on DIN EN ISO 7726.

    Based on the air velocity, this function will decide whether to use the
    natural or forced convection.

    Calculate hcg (the coefficient of heat transfer) for both natural and forced
    convection. Check which one is higher and use that (defined in Section B.2.3)

    This function performs better for larger arrays. For smaller arrays, the
    numpy-based function outperforms this function.

    :param ta: air temperature in °C
    :param tg: black globe temperature in °C
    :param v: air velocity in m/s
    :param d: diameter of the black globe in m (default 0.15 m)
    :param e: emissivity of the black globe (default 0.95)

    :returns: returns the mean radiant temperature in °C
    """
    tg = np.array(tg)
    v = np.array(v)
    ta = np.array(ta)
    # check if we're using the default values for d and e
    if isinstance(d, float):
        d = np.full_like(tg, d, dtype=float)
    else:
        d = np.array(d)
    if isinstance(e, float):
        e = np.full_like(tg, e, dtype=float)
    else:
        e = np.array(e)

    # 1. check for correct shape
    if not (
            tg.ndim <= 1 and v.ndim <= 1 and ta.ndim <= 1 and
            d.ndim <= 1 and e.ndim <= 1
    ):
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for same length
    if not (tg.size == v.size == ta.size == d.size == e.size):
        raise ValueError('All arrays must have the same length')

    # 3. check for value ranges
    if np.any(d <= 0):
        raise ValueError('The globe diameter (d) must be positive')

    if np.any((e < 0) | (e > 1)):
        raise ValueError('The emissivity (e) must be between 0 and 1')

    result = thermal_comfort_mod.mean_radiant_temp(tg=tg, v=v, ta=ta, d=d, e=e)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


# autopep8: off
@overload
def wet_bulb_temp(
        ta: float,
        rh: float,
) -> float: ...


@overload
def wet_bulb_temp(
        ta: Iterable[float],
        rh: Iterable[float],
) -> npt.NDArray[T]: ...


@overload
def wet_bulb_temp(
        ta: npt.NDArray[T],
        rh: npt.NDArray[T],
) -> npt.NDArray[T]: ...
# autopep8: on


def wet_bulb_temp(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
        rh: Union[npt.NDArray[T], float, Iterable[float]],
) -> Union[npt.NDArray[T], float]:
    """Calculate the wet bulb temperature following the Stull (2011) equation

    :param ta: air temperature in °C
    :param rh: relative humidity in %

    :returns: returns the wet bulb temperature in °C

    **References**

    - Stull, R., 2011. Wet-Bulb Temperature from Relative Humidity and
      Air Temperature. J. Appl. Meteorol. Climatol. 50, 2267-2269.
      https://doi.org/10.1175/JAMC-D-11-0143.1
    """
    ta = np.array(ta)
    rh = np.array(rh)

    # 1. check for correct shape
    if not (ta.ndim <= 1 and rh.ndim <= 1):
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for same length
    if not (ta.size == rh.size):
        raise ValueError('All arrays must have the same length')

    result = thermal_comfort_mod.wet_bulb_temp(ta=ta, rh=rh)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


# autopep8: off
@overload
def sat_vap_press_water(ta: float) -> float: ...


@overload
def sat_vap_press_water(ta: Iterable[float]) -> npt.NDArray[T]: ...


@overload
def sat_vap_press_water(ta: npt.NDArray[T]) -> npt.NDArray[T]: ...
# autopep8: on


def sat_vap_press_water(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
) -> Union[npt.NDArray[T], float]:
    """Calculate the saturation vapor pressure **over water** following the
    equation in VDI 3786 sheet 04.

    This functions is optimized on 1D-array operations, however also scalars may
    be provided.

    :param ta: air temperature in °C

    :returns: returns the saturation vapor pressure in hPa

    **References**

    - Sonntag, D. (1990). Important new values of the physical con- stants of 1986,
      vapour pressure formulations based on the ITC-90, and psychrometer formulae.
      Z. Meteorol., 40, 340–344.
    """
    ta = np.array(ta)

    # 1. check for correct shape
    if not ta.ndim <= 1:
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )

    # 2. check for value ranges
    if np.any((ta < -45) | (ta > 60)):
        raise ValueError('The air temperature (ta) must be between -45 and 60 °C')

    result = thermal_comfort_mod.sat_vap_press_water(ta=ta)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


# autopep8: off
@overload
def sat_vap_press_ice(ta: float) -> float: ...


@overload
def sat_vap_press_ice(ta: Iterable[float]) -> npt.NDArray[T]: ...


@overload
def sat_vap_press_ice(ta: npt.NDArray[T]) -> npt.NDArray[T]: ...
# autopep8: on


def sat_vap_press_ice(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
) -> Union[npt.NDArray[T], float]:
    """Calculate the saturation vapor pressure **over ice** following the
    equation in VDI 3786 sheet 04.

    This functions is optimized on 1D-array operations, however also scalars may
    be provided.

    :param ta: air temperature in °C

    :returns: returns the saturation vapor pressure in hPa

    **References**

    - Sonntag, D. (1990). Important new values of the physical con- stants of 1986,
      vapour pressure formulations based on the ITC-90, and psychrometer formulae.
      Z. Meteorol., 40, 340–344.
    """

    ta = np.array(ta)

    # 1. check for correct shape
    if not ta.ndim <= 1:
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for value ranges
    if np.any((ta < -65) | (ta > 0.01)):
        raise ValueError('The air temperature (ta) must be between -65 and 0.01 °C')

    result = thermal_comfort_mod.sat_vap_press_ice(ta=ta)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


# autopep8: off
@overload
def dew_point(
        ta: float,
        rh: float,
) -> float: ...


@overload
def dew_point(
        ta: Iterable[float],
        rh: Iterable[float],
) -> npt.NDArray[T]: ...


@overload
def dew_point(
        ta: npt.NDArray[T],
        rh: npt.NDArray[T],
) -> npt.NDArray[T]: ...
# autopep8: on


def dew_point(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
        rh: Union[npt.NDArray[T], float, Iterable[float]],
) -> Union[npt.NDArray[T], float]:
    """Calculate the dew point following the equation in VDI 3786 sheet 04.

    For temperatures values grate or equal to 0 °C, the dew point is calculated
    using the saturation vapor pressure over water. For temperatures below 0 °C,
    the dew point is calculated using the saturation vapor pressure over ice.

    This functions is optimized on 1D-array operations, however also scalars may
    be provided.

    :param ta: air temperature in °C
    :param rh: relative humidity in %

    :returns: returns the dew point temperature in °C

    **References**

    - Sonntag, D. (1990). Important new values of the physical con- stants of 1986,
      vapour pressure formulations based on the ITC-90, and psychrometer formulae.
      Z. Meteorol., 40, 340–344.
    - Sonntag, D. (1994). Advancements in the field of hygrometry.
      Meteorologische Zeitschrift, 51–66. https://doi.org/10.1127/metz/3/1994/51
    """
    ta = np.array(ta)
    rh = np.array(rh)

    # 1. check for correct shape
    if not (ta.ndim <= 1 and rh.ndim <= 1):
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for same length
    if not (ta.size == rh.size):
        raise ValueError('All arrays must have the same length')

    result = thermal_comfort_mod.dew_point(ta=ta, rh=rh)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


# autopep8: off
@overload
def absolute_humidity(
        ta: float,
        rh: float,
) -> float: ...


@overload
def absolute_humidity(
        ta: Iterable[float],
        rh: Iterable[float],
) -> npt.NDArray[T]: ...


@overload
def absolute_humidity(
        ta: npt.NDArray[T],
        rh: npt.NDArray[T],
) -> npt.NDArray[T]: ...
# autopep8: on


def absolute_humidity(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
        rh: Union[npt.NDArray[T], float, Iterable[float]],
) -> Union[npt.NDArray[T], float]:
    """Calculate the absolute humidity above water following the equation in
    VDI 3786 sheet 04.

    This functions is optimized on 1D-array operations, however also scalars may
    be provided.

    :param ta: air temperature in °C
    :param rh: relative humidity in %

    :returns: returns the absolute humidity in g/m³
    """
    ta = np.array(ta)
    rh = np.array(rh)

    # 1. check for correct shape
    if not (ta.ndim <= 1 and rh.ndim <= 1):
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for same length
    if not (ta.size == rh.size):
        raise ValueError('All arrays must have the same length')

    result = thermal_comfort_mod.absolute_humidity(ta=ta, rh=rh)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


# autopep8: off
@overload
def specific_humidity(
        ta: float,
        rh: float,
        p: float = 1013.25,
) -> float: ...


@overload
def specific_humidity(
        ta: Iterable[float],
        rh: Iterable[float],
        p: Union[Iterable[float], float] = 1013.25,
) -> npt.NDArray[T]: ...


@overload
def specific_humidity(
        ta: npt.NDArray[T],
        rh: npt.NDArray[T],
        p: Union[npt.NDArray[T], float] = 1013.25,
) -> npt.NDArray[T]: ...
# autopep8: on


def specific_humidity(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
        rh: Union[npt.NDArray[T], float, Iterable[float]],
        p: Union[npt.NDArray[T], float, Iterable[float]] = 1013.25,
) -> Union[npt.NDArray[T], float]:
    """Calculate the specific humidity above water following the equation in
    VDI 3786 sheet 04.

    This functions is optimized on 1D-array operations, however also scalars may
    be provided.

    :param ta: air temperature in °C
    :param rh: relative humidity in %
    :param p: air pressure in hPa (default 1013.25 hPa)

    :returns: returns the absolute humidity in g/kg
    """
    ta = np.array(ta)
    rh = np.array(rh)
    # when we use the default, we need to reshape the arrays
    if isinstance(p, float):
        p = np.full_like(ta, p, dtype=float)
    else:
        p = np.array(p)

    # 1. check for correct shape
    if not (ta.ndim <= 1 and rh.ndim <= 1 and p.ndim <= 1):
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for same length
    if not (ta.size == rh.size == p.size):
        raise ValueError('All arrays must have the same length')

    result = thermal_comfort_mod.specific_humidity(ta=ta, rh=rh, p=p)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result
