from typing import TypeVar

import numpy as np
import numpy.typing as npt

T = TypeVar('T', bound=np.floating | np.integer)


def _tmrt_forced_convection(
        tg: npt.NDArray[T] | float,
        ta: npt.NDArray[T] | float,
        va: npt.NDArray[T] | float,
        d: npt.NDArray[T] | float,
        e: npt.NDArray[T] | float,
) -> npt.NDArray[T] | np.floating:
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
        tg: npt.NDArray[T] | float,
        ta: npt.NDArray[T] | float,
        d: npt.NDArray[T] | float,
        e: npt.NDArray[T] | float,
) -> npt.NDArray[T] | np.floating:
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
        tg: npt.NDArray[T] | float,
        va: npt.NDArray[T] | float,
        ta: npt.NDArray[T] | float,
        d: npt.NDArray[T] | float = 0.15,
        e: npt.NDArray[T] | float = 0.95,
) -> npt.NDArray[T] | np.floating:
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
