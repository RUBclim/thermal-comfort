from typing import Any
from typing import TypeVar
from typing import Union

import numpy as np
import numpy.typing as npt

from ._thermal_comfort import thermal_comfort_mod


T = TypeVar('T', bound=Union[np.floating, np.integer])


def _tmrt_forced_convection(
        tg: Union[npt.NDArray[T], float],
        ta: Union[npt.NDArray[T], float],
        va: Union[npt.NDArray[T], float],
        d: Union[npt.NDArray[T], float],
        e: Union[npt.NDArray[T], float],
) -> Union[npt.NDArray[T], np.floating]:
    """
    Calculate the mean radiant temperature based on DIN EN ISO 7726 for forced
    convection.

    tg = black globe temperature
    ta = air temperature
    va = air velocity
    d = diameter of the black globe
    e = emissivity of the black globe
    """
    return (
        np.power(
            (
                (np.power((tg + 273), 4)) +
                (
                    ((1.1 * np.power(10, 8)) * np.power(va, 0.6)) /
                    (e * np.power(d, 0.4))
                ) * (tg - ta)
            ), 0.25,
        )
    ) - 273


def _tmrt_natural_convection(
        tg: Union[npt.NDArray[T], float],
        ta: Union[npt.NDArray[T], float],
        d: Union[npt.NDArray[T], float],
        e: Union[npt.NDArray[T], float],
) -> Union[npt.NDArray[T], np.floating]:
    """
    Calculate the mean radiant temperature based on DIN EN ISO 7726 for natural
    convection.

    tg = black globe temperature
    ta = air temperature
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
        tg: Union[npt.NDArray[T], float],
        va: Union[npt.NDArray[T], float],
        ta: Union[npt.NDArray[T], float],
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

    :param tg: black globe temperature
    :param va: air velocity
    :param ta: air temperature
    :param d: diameter of the black globe (default 0.15 m)
    :param e: emissivity of the black globe (default 0.95)
    """
    tg = np.asarray(tg)
    va = np.asarray(va)
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
    hcg_forced = (6.3 * (np.power(va, 0.6) / np.power(d, 0.4)))

    mask = hcg_natural > hcg_forced
    output = np.zeros(shape=mask.shape)
    # apply the natural convection cases
    output[mask] = _tmrt_natural_convection(
        tg=tg[mask], ta=ta[mask], d=d[mask], e=e[mask],
    )
    # apply the forced convection cases
    output[~mask] = _tmrt_forced_convection(
        tg=tg[~mask], ta=ta[~mask], va=va[~mask], d=d[~mask], e=e[~mask],
    )
    return output


def mrt(
        tg: npt.ArrayLike,
        va: npt.ArrayLike,
        ta: npt.ArrayLike,
        d: npt.ArrayLike = 0.15,
        e: npt.ArrayLike = 0.95,
) -> npt.NDArray[Any]:
    """
    Calculate the mean radiant temperature based on DIN EN ISO 7726.

    Based on the air velocity, this function will decide whether to use the
    natural or forced convection.

    Calculate hcg (the coefficient of heat transfer) for both natural and forced
    convection. Check which one is higher and use that (defined in Section B.2.3)

    This function performs better for larger arrays. For smaller arrays, the
    numpy-based function outperforms this function.

    :param tg: black globe temperature
    :param va: air velocity
    :param ta: air temperature
    :param d: diameter of the black globe (default 0.15 m)
    :param e: emissivity of the black globe (default 0.95)
    """
    tg = np.array(tg)
    va = np.array(va)
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
            tg.ndim <= 1 and va.ndim <= 1 and ta.ndim <= 1 and
            d.ndim <= 1 and e.ndim <= 1
    ):
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for same length
    if not (tg.size == va.size == ta.size == d.size == e.size):
        raise ValueError('All arrays must have the same length')

    # 3. check for value ranges
    if np.any(d <= 0):
        raise ValueError('The globe diameter (d) must be positive')

    if np.any((e < 0) | (e > 1)):
        raise ValueError('The emissivity (e) must be between 0 and 1')

    result = thermal_comfort_mod.mrt(tg=tg, va=va, ta=ta, d=d, e=e)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


def twb(
        ta: npt.ArrayLike,
        rh: npt.ArrayLike,
) -> npt.NDArray[Any]:
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
