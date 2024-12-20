from numpy.typing import ArrayLike
from typing import Any
import numpy.typing as npt


def pet_static(ta: float, rh: float, v:float, tmrt: float, p: float) -> float:
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

    The procedure does not work on arrays. It has to be wrappen in ``np.vectorize`` to
    properly work on arrays.

    :param ta: air temperature in °C
    :param rh: relative humidity in %
    :param v: wind speed in m/s
    :param tmrt: mean radiant temperature in °C
    :param p: atmospheric pressure in hPa

    :returns: Physiological Equivalent Temperature (PET) in °C
    """
    ...


def utci_approx(
        ta: ArrayLike,
        tmrt: ArrayLike,
        va:ArrayLike,
        rh:ArrayLike,
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
    :param va: wind speed in m/s
    :param rh: relative humidity in %

    :returns: Universal Thermal Climate Index (UTCI) in °C
    """
    ...
