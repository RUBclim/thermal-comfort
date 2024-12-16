from numpy.typing import ArrayLike
from typing import Any
import numpy.typing as npt

def utci_approx(ta: float, tmrt: float, va:float, rh:float) -> float:
    """test 123"""
    ...

def es(ta: float) -> float:
    """Calculates saturation vapour pressure over water in hPa for input air
    temperature (ta) in celsius according to:
    Hardy, R.; ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature,
    Dewpoint Temperature and Enhancement Factors in the Range -100 to 100 °C;
    Proceedings of Third International Symposium on Humidity and Moisture;
    edited by National Physical Laboratory (NPL), London, 1998, pp. 214-221
    http://www.thunderscientific.com/tech_info/reflibrary/its90formulas.pdf
    (retrieved 2008-10-01)

    :param ta: air temperature in °C

    :return: saturation vapour pressure over water in hPa
    """
    ...


def es_vectorized(ta: ArrayLike) -> npt.NDArray[Any]: ...


def utci_approx_vectorized(ta: ArrayLike, tmrt: ArrayLike, va:ArrayLike, rh:ArrayLike) -> npt.NDArray[Any]:
    """This will flatten n-dimensional arrays to 1D arrays"""
    ...
