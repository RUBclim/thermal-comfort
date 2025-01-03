import warnings
from typing import Any

import numpy as np
import numpy.typing as npt

from ._thermal_comfort import thermal_comfort_mod


def utci_approx(
        ta: npt.ArrayLike,
        tmrt: npt.ArrayLike,
        va: npt.ArrayLike,
        rh: npt.ArrayLike,
) -> npt.NDArray[Any]:
    """Calculate the Universal Thermal Climate Index (UTCI)

    The UTCI is implemented as described in VDI 3787 Part 2. The fortran code was
    vendored from here:

    - https://utci.org/resources/UTCI%20Program%20Code.zip


    The procedure works on 1D-arrays.

    :param ta: air temperature in °C
    :param tmrt: mean radiant temperature in °C
    :param va: wind speed in m/s
    :param rh: relative humidity in %

    :returns: Universal Thermal Climate Index (UTCI) in °C
    """
    ta = np.array(ta)
    tmrt = np.array(tmrt)
    va = np.array(va)
    rh = np.array(rh)

    # 1. check for correct shape
    if not (ta.ndim == 1 and tmrt.ndim == 1 and va.ndim == 1 and rh.ndim == 1):
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for same length
    if not (ta.size == tmrt.size == va.size == rh.size):
        raise ValueError('All arrays must have the same length')

    # 3. check for value ranges
    if np.any((ta < -50) | (ta > 50)):
        warnings.warn(
            'encountered a value for ta outside of the defined range of '
            '-50 <= ta <= 50 °C',
            category=RuntimeWarning,
            stacklevel=2,
        )
    if np.any((va < 0.5) | (va > 17)):
        warnings.warn(
            'encountered a value for va outside of the defined range of '
            '0.5 <= va <= 17',
            stacklevel=2,
            category=RuntimeWarning,
        )
    delta_ta_tmrt = tmrt - ta
    if np.any((delta_ta_tmrt < -30) | (delta_ta_tmrt > 70)):
        warnings.warn(
            'encountered a value for tmrt outside of the defined range of '
            '-30 °C below or 70 °C above ta',
            category=RuntimeWarning,
            stacklevel=2,
        )

    return thermal_comfort_mod.utci_approx(ta=ta, tmrt=tmrt, va=va, rh=rh)


def pet_static(
        ta: npt.ArrayLike,
        rh: npt.ArrayLike,
        v: npt.ArrayLike,
        tmrt: npt.ArrayLike,
        p: npt.ArrayLike,
) -> npt.NDArray[Any]:
    """Calculate the Physiological Equivalent Temperature (PET).

    The PET is implemented as described in VDI 3787 Part 2. The fortran code was
    vendored from here:

    - https://www.vdi.de/richtlinien/programme-zu-vdi-richtlinien/vdi-3787-blatt-2
    - http://web.archive.org/web/20241219155627/https://www.vdi.de/fileadmin/pages/vdi_de/redakteure/ueber_uns/fachgesellschaften/KRdL/dateien/VDI_3787-2_PET.zip

    The code was adapted to retrieve relative humidity instead of vapor pressure. The
    saturation vapor pressure is calculated using the Wexler formula.

    The procedure has some limitations compare to a full implementation of the PET.
    Many variables are set to static values, such as:

    - ``age = 35.``
    - ``mbody = 75.``
    - ``ht = 1.75``
    - ``work = 80.``
    - ``eta = 0.``
    - ``icl = 0.9``
    - ``fcl = 1.15``
    - ``pos = 1``
    - ``sex = 1``

    :param ta: air temperature in °C
    :param rh: relative humidity in %
    :param v: wind speed in m/s
    :param tmrt: mean radiant temperature in °C
    :param p: atmospheric pressure in hPa

    :returns: Physiological Equivalent Temperature (PET) in °C
    """  # noqa: E501
    ta = np.array(ta)
    rh = np.array(rh)
    v = np.array(v)
    tmrt = np.array(tmrt)
    p = np.array(p)

    # 1. check for correct shape
    if not (
            ta.ndim == 1 and rh.ndim == 1 and v.ndim == 1 and tmrt.ndim == 1
            and p.ndim == 1
    ):
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for same length
    if not (ta.size == rh.size == v.size == tmrt.size):
        raise ValueError('All arrays must have the same length')

    return thermal_comfort_mod.pet_static(ta=ta, rh=rh, v=v, tmrt=tmrt, p=p)
