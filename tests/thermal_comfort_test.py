import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from thermal_comfort import pet_static
from thermal_comfort import utci_approx


def load_test_data():
    with open('testing/utci_references.txt') as f:
        # skip the headers
        lines = f.readlines()[35:]
        values = []
        for line in lines:
            values.append([float(i) for i in line.split('\t')])

    return values


@pytest.mark.parametrize(
    (
        'ta', 'd_tmrt', 'va', 'rh', 'pa', 'offset',
        'utci', 'utci_table', 'utci_polynomial',
    ),
    (*load_test_data(),),
)
def test_utci_approx(
        ta,
        d_tmrt,
        va,
        rh,
        pa,
        offset,
        utci,
        utci_table,
        utci_polynomial,
):
    assert pytest.approx(
        utci_approx(ta=ta, tmrt=d_tmrt + ta, va=va, rh=rh),
        abs=1e-1,
    ) == utci_polynomial


def test_utci_approx_native_vectorized():
    data = np.array(load_test_data())
    ta = data[:, 0]
    tmrt = data[:, 1] + data[:, 0]
    va = data[:, 2]
    rh = data[:, 3]
    expected = data[:, 8]

    assert_array_almost_equal(
        utci_approx(ta=ta, tmrt=tmrt, va=va, rh=rh),
        expected,
        decimal=1,
    )


def test_utci_approx_native_vectorized_2d_array():
    ta = np.array([[-49.9, -49.8], [-49.2, -49]])
    tmrt = np.array([[-16.7 + -49.9, -14.6 + -49.8], [-5.6 + -49.2, -17.2 + -49]])
    va = np.array([[28, 14], [18, 5]])
    rh = np.array([[78, 77], [100, 98]])
    expected = np.array([[-95.2, -90.4], [-98.7, -66.1]])

    assert_array_almost_equal(
        utci_approx(
            ta=ta,
            tmrt=tmrt,
            va=va,
            rh=rh,
        ).reshape((2, 2), order='F'),
        expected,
        decimal=1,
    )


def _rh(vpa, ta):
    """we changed the interface of pet to use relative humidity instead of
    vapour pressure, hence we need to calculate the relative humidity from
    the vapour pressure and air temperature since only this is specified in the
    paper by Höppe 1999
    """
    return (vpa * 100) / (6.1094 * np.exp(17.625 * ta / (ta + 243.04)))


@pytest.mark.parametrize(
    ('ta', 'rh', 'v', 'tmrt', 'expected'),
    # reference values are from Höppe 1999
    (
        (21, _rh(12, 21), 0.1, 21, 21),
        (-5, _rh(2, -5), 0.5, 40, 10),
        (-5, _rh(2, -5), 6, -5, -13),
        (30, _rh(21, 30), 1, 60, 43),
        (30, _rh(21, 30), 1, 30, 29),
    ),
)
def test_pet_static(ta, rh, v, tmrt, expected):
    # the values are only supplied as integers
    assert np.round(
        pet_static(ta=ta, rh=rh, v=v, tmrt=tmrt, p=1013.25),
    ) == expected
