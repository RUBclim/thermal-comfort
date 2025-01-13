from typing import overload
from typing import TypeVar
from typing import Union

import numpy as np
import numpy.typing as npt

from ._thermal_comfort import thermal_comfort_mod


T = TypeVar('T', bound=Union[np.floating, np.integer])


def _tmrt_forced_convection(
        ta: Union[npt.NDArray[T], float],
        tg: Union[npt.NDArray[T], float],
        v: Union[npt.NDArray[T], float],
        d: Union[npt.NDArray[T], float],
        e: Union[npt.NDArray[T], float],
) -> Union[npt.NDArray[T], np.floating]:
    """
    Calculate the mean radiant temperature based on DIN EN ISO 7726 for forced
    convection.

    ta = air temperature
    tg = black globe temperature
    v = air velocity
    d = diameter of the black globe
    e = emissivity of the black globe
    """
    return (
        np.power(
            (
                (np.power((tg + 273), 4)) +
                (
                    ((1.1 * np.power(10, 8)) * np.power(v, 0.6)) /
                    (e * np.power(d, 0.4))
                ) * (tg - ta)
            ), 0.25,
        )
    ) - 273


def _tmrt_natural_convection(
        ta: Union[npt.NDArray[T], float],
        tg: Union[npt.NDArray[T], float],
        d: Union[npt.NDArray[T], float],
        e: Union[npt.NDArray[T], float],
) -> Union[npt.NDArray[T], np.floating]:
    """
    Calculate the mean radiant temperature based on DIN EN ISO 7726 for natural
    convection.

    ta = air temperature
    tg = black globe temperature
    d = diameter of the black globe
    e = emissivity of the black globe
    """
    return (
        np.power(
            ((np.power(tg + 273, 4)) + (np.power(0.25*10, 8)/e) *
             (np.power((np.abs(tg - ta)/d), 0.25))*(tg - ta)),
            0.25,
        )
    ) - 273


def mrt_np(
        ta: Union[npt.NDArray[T], float],
        tg: Union[npt.NDArray[T], float],
        v: Union[npt.NDArray[T], float],
        d: Union[npt.NDArray[T], float] = 0.15,
        e: Union[npt.NDArray[T], float] = 0.95,
) -> Union[npt.NDArray[T], np.floating]:
    """
    Calculate the mean radiant temperature based on DIN EN ISO 7726.

    Based on the air velocity, this function will decide whether to use the
    natural or forced convection.

    Calculate hcg (the coefficient of heat transfer) for both natural and forced
    convection. Check which one is higher and use that (defined in Section B.2.3)

    This function performs better for smaller arrays. For larger arrays, the
    fortran-based function outperforms this function.

    :param ta: air temperature in °C
    :param tg: black globe temperature in °C
    :param v: air velocity in m/s
    :param d: diameter of the black globe in m (default 0.15 m)
    :param e: emissivity of the black globe (default 0.95)
    """
    tg = np.asarray(tg)
    v = np.asarray(v)
    ta = np.asarray(ta)
    # when we use the default, we need to reshape the arrays
    if isinstance(d, float):
        d = np.full_like(tg, d, dtype=float)
    if isinstance(e, float):
        e = np.full_like(tg, e, dtype=float)

    # check for value ranges
    if np.any(d <= 0):
        raise ValueError('The globe diameter (d) must be positive')

    if np.any((e < 0) | (e > 1)):
        raise ValueError('The emissivity (e) must be between 0 and 1')

    # we need to mask the values and apply the correct function to these values
    hcg_natural = np.power(1.4 * (np.abs(tg - ta) / 0.15), 0.25)
    hcg_forced = (6.3 * (np.power(v, 0.6) / np.power(d, 0.4)))

    mask = hcg_natural > hcg_forced
    output = np.zeros(shape=mask.shape)
    # apply the natural convection cases
    output[mask] = _tmrt_natural_convection(
        tg=tg[mask], ta=ta[mask], d=d[mask], e=e[mask],
    )
    # apply the forced convection cases
    output[~mask] = _tmrt_forced_convection(
        tg=tg[~mask], ta=ta[~mask], v=v[~mask], d=d[~mask], e=e[~mask],
    )
    return output


# autopep8: off
@overload
def mrt(
        ta: float,
        tg: float,
        v: float,
        d: float = 0.15,
        e: float = 0.95,
) -> float: ...


@overload
def mrt(
        ta: npt.NDArray[T],
        tg: npt.NDArray[T],
        v: npt.NDArray[T],
        d: Union[npt.NDArray[T], float] = 0.15,
        e: Union[npt.NDArray[T], float] = 0.95,
) -> npt.NDArray[T]: ...
# autopep8: on


def mrt(
        ta: Union[npt.NDArray[T], float],
        tg: Union[npt.NDArray[T], float],
        v: Union[npt.NDArray[T], float],
        d: Union[npt.NDArray[T], float] = 0.15,
        e: Union[npt.NDArray[T], float] = 0.95,
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

    result = thermal_comfort_mod.mrt(tg=tg, v=v, ta=ta, d=d, e=e)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


# autopep8: off
@overload
def twb(
        ta: float,
        rh: float,
) -> float: ...


@overload
def twb(
        ta: npt.NDArray[T],
        rh: npt.NDArray[T],
) -> npt.NDArray[T]: ...
# autopep8: on


def twb(
        ta: Union[npt.NDArray[T], float],
        rh: Union[npt.NDArray[T], float],
) -> Union[npt.NDArray[T], float]:
    """Calculate the wet bulb temperature following the Stull (2011) equation

    :param ta: air temperature in °C
    :param rh: relative humidity in %

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

    result = thermal_comfort_mod.twb(ta=ta, rh=rh)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result
