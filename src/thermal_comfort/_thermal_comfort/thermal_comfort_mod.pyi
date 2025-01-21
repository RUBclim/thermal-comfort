from typing import Any
from typing import TypeVar

import numpy as np
import numpy.typing as npt

T = TypeVar('T', bound=np.floating | np.integer)


def pet_static(
        ta: npt.NDArray[T] | float,
        rh: npt.NDArray[T] | float,
        v: npt.NDArray[T] | float,
        tmrt: npt.NDArray[T] | float,
        p: npt.NDArray[T] | float,
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
    ...


def utci_approx(
        ta: npt.NDArray[T] | float,
        tmrt: npt.NDArray[T] | float,
        v: npt.NDArray[T] | float,
        rh: npt.NDArray[T] | float,
) -> npt.NDArray[Any]:
    """Calculate the Universal Thermal Climate Index (UTCI)

    The UTCI is implemented as described in VDI 3787 Part 2. The fortran code was
    vendored from here:

    - https://utci.org/resources/UTCI%20Program%20Code.zip


    The procedure works on 1D-arrays. Higher dimensional arrays are flattened into
    1D-arrays. The output array is always 1D. You can reshape the output array back to
    its original shape. using: ``.reshape(<shape>, order='F')``. This function uses
    vectors to improve performance.

    :param ta: air temperature in °C
    :param tmrt: mean radiant temperature in °C
    :param v: wind speed in m/s
    :param rh: relative humidity in %

    :returns: Universal Thermal Climate Index (UTCI) in °C
    """
    ...


def mean_radiant_temp(
        ta: npt.NDArray[T] | float,
        tg: npt.NDArray[T] | float,
        v: npt.NDArray[T] | float,
        d: npt.NDArray[T] | float,
        e: npt.NDArray[T] | float,
) -> npt.NDArray[T]:
    """
    Calculate the mean radiant temperature based on DIN EN ISO 7726.

    Based on the air velocity, this function will decide whether to use the
    natural or forced convection.

    Calculate hcg (the coefficient of heat transfer) for both natural and forced
    convection. Check which one is higher and use that (defined in Section B.2.3)

    This function performs better for larger arrays. For smaller arrays, the
    numpy-based function outperforms this function.

    :param ta: air temperature
    :param tg: black globe temperature
    :param v: air velocity
    :param d: diameter of the black globe (default 0.15 m)
    :param e: emissivity of the black globe (default 0.95)
    """
    ...


def wet_bulb_temp(
        ta: npt.NDArray[T] | float,
        rh: npt.NDArray[T] | float,
) -> npt.NDArray[T]:
    ...


def heat_index(
        ta: npt.NDArray[T] | float,
        rh: npt.NDArray[T] | float,
) -> npt.NDArray[T]:
    ...


def heat_index_extended(
        ta: npt.NDArray[T] | float,
        rh: npt.NDArray[T] | float,
) -> npt.NDArray[T]:
    ...


def sat_vap_press_water(ta: npt.NDArray[T] | float) -> npt.NDArray[T]: ...


def sat_vap_press_ice(ta: npt.NDArray[T] | float) -> npt.NDArray[T]: ...


def dew_point(
        ta: npt.NDArray[T] | float,
        rh: npt.NDArray[T] | float,
) -> npt.NDArray[T]:
    ...


def absolute_humidity(
        ta: npt.NDArray[T] | float,
        rh: npt.NDArray[T] | float,
) -> npt.NDArray[T]:
    ...


def specific_humidity(
        ta: npt.NDArray[T] | float,
        rh: npt.NDArray[T] | float,
        p: npt.NDArray[T] | float,
) -> npt.NDArray[T]:
    ...
