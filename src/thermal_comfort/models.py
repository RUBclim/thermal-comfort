import warnings
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
def utci_approx(
        ta: float,
        tmrt: float,
        v: float,
        rh: float,
) -> float: ...


@overload
def utci_approx(
        ta: Iterable[float],
        tmrt: Iterable[float],
        v: Iterable[float],
        rh: Iterable[float],
) -> npt.NDArray[T]: ...


@overload
def utci_approx(
        ta: npt.NDArray[T],
        tmrt: npt.NDArray[T],
        v: npt.NDArray[T],
        rh: npt.NDArray[T],
) -> npt.NDArray[T]: ...
# autopep8: on


def utci_approx(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
        tmrt: Union[npt.NDArray[T], float, Iterable[float]],
        v: Union[npt.NDArray[T], float, Iterable[float]],
        rh: Union[npt.NDArray[T], float, Iterable[float]],
) -> Union[npt.NDArray[T], float]:
    """Calculate the Universal Thermal Climate Index (UTCI).

    The UTCI is implemented as described in VDI 3787 Part 2. The fortran code was
    vendored from here: https://utci.org/resources/UTCI%20Program%20Code.zip.
    A few changes were implemented to the original fortran-code.

    - Instead of taking the vapor pressure as an argument, it now takes the
      relative humidity. The vapor pressure is calculated from the relative humidity
      using the formula by Wexler (1976) which is described by Hardy (1998).
    - Arguments were renamed for a consistent interface within this package.


    This functions is optimized on 1D-array operations, however also scalars may
    be provided.

    :param ta: Air temperature in °C in the range between -50 and 50 °C
    :param tmrt: Mean radiant temperature in °C in the range between -30 °C below
        and 70 °C above ta
    :param v: Wind speed in m/s in the range between 0.5 and 17 m/s
    :param rh: Relative humidity in %

    :returns: Universal Thermal Climate Index (UTCI) in °C

    **References**

    - Hardy, R.; ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature,
      Dewpoint Temperature and Enhancement Factors in the Range -100 to 100 °C;
      Proceedings of Third International Symposium on Humidity and Moisture;
      edited by National Physical Laboratory (NPL), London, 1998, pp. 214-221
      https://www.decatur.de/javascript/dew/resources/its90formulas.pdf
    - Wexler, A., Vapor Pressure Formulation for Water in Range 0 to 100°C.
      A Revision, Journal of Research of the National Bureau of Standards - A.
      Physics and Chemistry, September - December 1976, Vol. 80A, Nos. 5 and 6, 775-78
    """
    ta = np.array(ta)
    tmrt = np.array(tmrt)
    v = np.array(v)
    rh = np.array(rh)

    # 1. check for correct shape
    if not (ta.ndim <= 1 and tmrt.ndim <= 1 and v.ndim <= 1 and rh.ndim <= 1):
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for same length
    if not (ta.size == tmrt.size == v.size == rh.size):
        raise ValueError('All arrays must have the same length')

    # 3. check for value ranges
    if np.any((ta < -50) | (ta > 50)):
        warnings.warn(
            'encountered a value for ta outside of the defined range of '
            '-50 <= ta <= 50 °C',
            category=RuntimeWarning,
            stacklevel=2,
        )
    if np.any((v < 0.5) | (v > 17)):
        warnings.warn(
            'encountered a value for v outside of the defined range of '
            '0.5 <= v <= 17',
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

    result = thermal_comfort_mod.utci_approx(ta=ta, tmrt=tmrt, v=v, rh=rh)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


# autopep8: off
@overload
def pet_static(
        ta: float,
        tmrt: float,
        v: float,
        rh: float,
        p: float,
) -> float: ...


@overload
def pet_static(
        ta: Iterable[float],
        tmrt: Iterable[float],
        v: Iterable[float],
        rh: Iterable[float],
        p: Iterable[float],
) -> npt.NDArray[T]: ...


@overload
def pet_static(
        ta: npt.NDArray[T],
        tmrt: npt.NDArray[T],
        v: npt.NDArray[T],
        rh: npt.NDArray[T],
        p: npt.NDArray[T],
) -> npt.NDArray[T]: ...
# autopep8: on


def pet_static(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
        tmrt: Union[npt.NDArray[T], float, Iterable[float]],
        v: Union[npt.NDArray[T], float, Iterable[float]],
        rh: Union[npt.NDArray[T], float, Iterable[float]],
        p: Union[npt.NDArray[T], float, Iterable[float]],
) -> Union[npt.NDArray[T], float]:
    """Calculate the Physiological Equivalent Temperature (PET).

    The PET is implemented as described in VDI 3787 Part 2. The fortran code was
    vendored from here:

    - https://www.vdi.de/richtlinien/programme-zu-vdi-richtlinien/vdi-3787-blatt-2
    - http://web.archive.org/web/20241219155627/https://www.vdi.de/fileadmin/pages/vdi_de/redakteure/ueber_uns/fachgesellschaften/KRdL/dateien/VDI_3787-2_PET.zip

    - Instead of taking the vapor pressure as an argument, it now takes the
      relative humidtiy. The vapor pressure is calculated from the relative humidtiy
      using the formular by Wexler (1976) which is described by Hardy (1998).
    - Arguments were renamed for a consistent interface within this package.

    This functions is optimized on 1D-array operations, however also scalars may be provided.

    The procedure has some limitations compared to a full implementation of the PET.
    Many variables are set to static values, such as:

    - ``age = 35``
    - ``mbody = 75``
    - ``ht = 1.75``
    - ``work = 80``
    - ``eta = 0``
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

    **References**

    - Hardy, R.; ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature, Dewpoint Temperature and Enhancement Factors in the Range -100 to 100 °C;
      Proceedings of Third International Symposium on Humidity and Moisture; edited by National Physical Laboratory (NPL), London, 1998, pp. 214-221
      https://www.decatur.de/javascript/dew/resources/its90formulas.pdf
    - Wexler, A., Vapor Pressure Formulation for Water in Range 0 to 100°C. A Revision, Journal of Research of
      the National Bureau of Standards - A. Physics and Chemistry, September - December 1976, Vol. 80A, Nos.
      5 and 6, 775-78
    """  # noqa: E501
    ta = np.array(ta)
    rh = np.array(rh)
    v = np.array(v)
    tmrt = np.array(tmrt)
    p = np.array(p)

    # 1. check for correct shape
    if not (
            ta.ndim <= 1 and rh.ndim <= 1 and v.ndim <= 1 and tmrt.ndim <= 1
            and p.ndim <= 1
    ):
        raise TypeError(
            'Only arrays with one dimension are allowed. '
            'Please reshape your array accordingly',
        )
    # 2. check for same length
    if not (ta.size == rh.size == v.size == tmrt.size == p.size):
        raise ValueError('All arrays must have the same length')

    result = thermal_comfort_mod.pet_static(ta=ta, rh=rh, v=v, tmrt=tmrt, p=p)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


# autopep8: off
@overload
def heat_index(ta: float, rh: float) -> float: ...


@overload
def heat_index(
        ta: Iterable[float],
        rh: Iterable[float],
) -> npt.NDArray[T]: ...


@overload
def heat_index(
        ta: npt.NDArray[T],
        rh: npt.NDArray[T],
) -> npt.NDArray[T]: ...
# autopep8: on


def heat_index(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
        rh: Union[npt.NDArray[T], float, Iterable[float]],
) -> Union[npt.NDArray[T], float]:
    r"""Calculate the heat index follwing Steadman R.G (1979) & Rothfusz L.P (1990).

    This calculates the heat index in the range of :math:`\ge` 80 °F (26.666 °C) and
    :math:`\ge` 40 % relative humidity. If values outside of this range are provided,
    they are returned as ``NaN``. This version natively works with °C as shown by
    Blazejczyk et al. 2011.

    This functions is optimized on 1D-array operations, however also scalars may
    be provided.

    :param ta: Air temperature in °C
    :param rh: Relative Humidity in %

    :returns: Heat Index in °C

    **References**

    - Steadman R.G. The assessment of sultriness. Part I: A temperature-humidity index
      based on human physiology and clothing. J. Appl. Meteorol. 1979;18:861-873.
      https://doi.org/10.1175/1520-0450(1979)018%3C0861:TAOSPI%3E2.0.CO;2

    - Rothfusz LP (1990) The heat index equation.
      NWS Southern Region Technical Attachment, SR/SSD 90-23, Fort Worth, Texas

    - Blazejczyk, K., Epstein, Y., Jendritzky, G. Staiger, H., Tinz, B. (2011)
      Comparison of UTCI to selected thermal indices. Int J Biometeorol 56,
      515-535 (2012). https://doi.org/10.1007/s00484-011-0453-2
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

    result = thermal_comfort_mod.heat_index(ta=ta, rh=rh)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result


# autopep8: off
@overload
def heat_index_extended(ta: float, rh: float) -> float: ...


@overload
def heat_index_extended(
        ta: Iterable[float],
        rh: Iterable[float],
) -> npt.NDArray[T]: ...


@overload
def heat_index_extended(
        ta: npt.NDArray[T],
        rh: npt.NDArray[T],
) -> npt.NDArray[T]: ...
# autopep8: on


def heat_index_extended(
        ta: Union[npt.NDArray[T], float, Iterable[float]],
        rh: Union[npt.NDArray[T], float, Iterable[float]],
) -> Union[npt.NDArray[T], float]:
    """Calculate the heat index following Steadman R.G (1979) & Rothfusz L.P (1990),
    but extends the range following The National Weather Service Weather Prediction
    Center.

    This function works for the entire range of air temperatures and relative humidity.
    It uses °F under the hood, since corrections are provided in °F only.

    This functions is optimized on 1D-array operations, however also scalars may
    be provided.

    :param ta: Air temperature in °C
    :param rh: Relative Humidity in %

    :returns: Heat Index in °C

    **References**

    - Steadman R.G. The assessment of sultriness. Part I: A temperature-humidity index
      based on human physiology and clothing. J. Appl. Meteorol. 1979;18:861-873.
      https://doi.org/10.1175/1520-0450(1979)018%3C0861:TAOSPI%3E2.0.CO;2

    - Rothfusz LP (1990) The heat index equation.
      NWS Southern Region Technical Attachment, SR/SSD 90-23, Fort Worth, Texas

    - NOAA/National Weather Service
      National Centers for Environmental Prediction Weather Prediction Center.
      The Heat Index Equation.
      https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml
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

    result = thermal_comfort_mod.heat_index_extended(ta=ta, rh=rh)
    # check if we have a single value
    if result.size == 1:
        return result.item()
    else:
        return result
