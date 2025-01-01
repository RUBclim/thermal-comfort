import math

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_array_almost_equal
from pandas.testing import assert_series_equal
from thermal_comfort import mrt
from thermal_comfort import mrt_np
from thermal_comfort import pet_static
from thermal_comfort import twb
from thermal_comfort import utci_approx


def load_test_data():
    with open('testing/utci_references.txt') as f:
        # skip the headers
        lines = f.readlines()[35:]
        values = []
        for line in lines:
            values.append([float(i) for i in line.split('\t')])

    return values


@pytest.mark.filterwarnings('ignore:encountered a value for')
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


@pytest.mark.parametrize(
    ('ta', 'tmrt', 'va', 'rh'),
    (
        (20, 20, 0.5, float('nan')),
        (20, 20, float('nan'), 50),
        (20, 20, float('nan'), float('nan')),
        (20, float('nan'), 0.5, 50),
        (20, float('nan'), 0.5, float('nan')),
        (20, float('nan'), float('nan'), 50),
        (20, float('nan'), float('nan'), float('nan')),
        (float('nan'), 20, 0.5, 50),
        (float('nan'), 20, 0.5, float('nan')),
        (float('nan'), 20, float('nan'), 50),
        (float('nan'), 20, float('nan'), float('nan')),
        (float('nan'), float('nan'), 0.5, 50),
        (float('nan'), float('nan'), 0.5, float('nan')),
        (float('nan'), float('nan'), float('nan'), 50),
        (float('nan'), float('nan'), float('nan'), float('nan')),
    ),
)
def test_utci_approx_missing_value(ta, tmrt, va, rh):
    assert math.isnan(utci_approx(ta=ta, tmrt=tmrt, va=va, rh=rh)[0])


@pytest.mark.filterwarnings('ignore:encountered a value for')
def test_utci_approx_with_vectors():
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


@pytest.mark.filterwarnings('ignore:encountered a value for')
def test_utci_approx_with_pandas_series():
    df = pd.DataFrame(
        load_test_data(),
        columns=[
            'Ta', 'Tr-Ta', 'va', 'rH', 'pa',
            'Offset', 'UTCI', 'UTCI_Table', 'UTCI_polynomial',
        ],
    )
    df['tmrt'] = df['Ta'] + df['Tr-Ta']
    df['utci_calc'] = utci_approx(
        ta=df['Ta'], tmrt=df['tmrt'], va=df['va'], rh=df['rH'],
    )
    assert_series_equal(
        left=df['UTCI_polynomial'],
        right=df['utci_calc'],
        atol=1e-1,
        check_names=False,
    )


def test_utci_approx_native_vectorized_2d_array():
    ta = np.array([[-49.9, -49.8], [-49.2, -49]])
    tmrt = np.array([[-16.7 + -49.9, -14.6 + -49.8], [-5.6 + -49.2, -17.2 + -49]])
    va = np.array([[8, 4], [8, 5]])
    rh = np.array([[78, 77], [100, 98]])
    expected = np.array([[-75.7, -63.5], [-73.9, -66.1]])

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


@pytest.mark.parametrize('ta', (-50.1, 50.1))
def test_utci_approx_warning_raised_when_outside_of_range_ta(ta):
    with pytest.warns(RuntimeWarning) as w:
        utci_approx(ta=ta, tmrt=60, va=1, rh=50)

    assert w[0].message.args[0] == (
        'encountered a value for ta outside of the defined range of -50 <= ta <= 50 °C'
    )


@pytest.mark.parametrize('tmrt', (-10.1, 90.1))
def test_utci_approx_warning_raised_when_outside_of_range_tmrt(tmrt):
    with pytest.warns(RuntimeWarning) as w:
        utci_approx(ta=20, tmrt=tmrt, va=1, rh=50)

    assert w[0].message.args[0] == (
        'encountered a value for tmrt outside of the defined range of '
        '-30 °C below or 70 °C above ta'
    )


@pytest.mark.parametrize('va', (0.49, 17.1))
def test_utci_approx_warning_raised_when_outside_of_range_va(va):
    with pytest.warns(RuntimeWarning) as w:
        utci_approx(ta=20, tmrt=60, va=va, rh=50)

    assert w[0].message.args[0] == (
        'encountered a value for va outside of the defined range of 0.5 <= va <= 17'
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


@pytest.mark.parametrize(
    ('ta', 'rh', 'v', 'tmrt', 'p'),
    (
        (20, 50, 0.2, 31, float('nan')),
        (20, 50, 0.2, float('nan'), 1013.25),
        (20, 50, 0.2, float('nan'), float('nan')),
        (20, 50, float('nan'), 31, 1013.25),
        (20, 50, float('nan'), 31, float('nan')),
        (20, 50, float('nan'), float('nan'), 1013.25),
        (20, 50, float('nan'), float('nan'), float('nan')),
        (20, float('nan'), 0.2, 31, 1013.25),
        (20, float('nan'), 0.2, 31, float('nan')),
        (20, float('nan'), 0.2, float('nan'), 1013.25),
        (20, float('nan'), 0.2, float('nan'), float('nan')),
        (20, float('nan'), float('nan'), 31, 1013.25),
        (20, float('nan'), float('nan'), 31, float('nan')),
        (20, float('nan'), float('nan'), float('nan'), 1013.25),
        (20, float('nan'), float('nan'), float('nan'), float('nan')),
        (float('nan'), 50, 0.2, 31, 1013.25),
        (float('nan'), 50, 0.2, 31, float('nan')),
        (float('nan'), 50, 0.2, float('nan'), 1013.25),
        (float('nan'), 50, 0.2, float('nan'), float('nan')),
        (float('nan'), 50, float('nan'), 31, 1013.25),
        (float('nan'), 50, float('nan'), 31, float('nan')),
        (float('nan'), 50, float('nan'), float('nan'), 1013.25),
        (float('nan'), 50, float('nan'), float('nan'), float('nan')),
        (float('nan'), float('nan'), 0.2, 31, 1013.25),
        (float('nan'), float('nan'), 0.2, 31, float('nan')),
        (float('nan'), float('nan'), 0.2, float('nan'), 1013.25),
        (float('nan'), float('nan'), 0.2, float('nan'), float('nan')),
        (float('nan'), float('nan'), float('nan'), 31, 1013.25),
        (float('nan'), float('nan'), float('nan'), 31, float('nan')),
        (float('nan'), float('nan'), float('nan'), float('nan'), 1013.25),
        (float('nan'), float('nan'), float('nan'), float('nan'), float('nan')),
    ),
)
def test_pet_static_missing_value(ta, rh, v, tmrt, p):
    assert math.isnan(pet_static(ta=ta, rh=rh, v=v, tmrt=tmrt, p=p))


def test_pet_static_missing_value_mixed_array():
    f = np.vectorize(pet_static, otypes=[float], cache=True)
    result = f(
        ta=np.array([20, float('nan')]),
        rh=np.array([50, float('nan')]),
        v=np.array([0.5, float('nan')]),
        tmrt=np.array([20, float('nan')]),
        p=np.array([1013.5, float('nan')]),
    )
    assert_array_almost_equal(result, [18, float('nan')], decimal=1)


def test_pet_static_missing_value_mixed_series():
    f = np.vectorize(pet_static, otypes=[float], cache=True)
    df = pd.DataFrame({
        'ta': pd.Series([20, float('nan')]),
        'rh': pd.Series([50, float('nan')]),
        'v': pd.Series([0.5, float('nan')]),
        'tmrt': pd.Series([20, float('nan')]),
        'p': pd.Series([1013.5, float('nan')]),
    })
    df['result'] = f(ta=df['ta'], rh=df['rh'], v=df['v'], tmrt=df['tmrt'], p=df['p'])
    assert_series_equal(
        left=df['result'],
        right=pd.Series([18, float('nan')]),
        atol=1,
        check_names=False,
    )


@pytest.mark.parametrize(
    ('ta', 'rh', 'v', 'tmrt', 'p'),
    (
        (20, 50, 0.2, 31, None),
        (20, 50, 0.2, None, 1013.25),
        (20, 50, 0.2, None, None),
        (20, 50, None, 31, 1013.25),
        (20, 50, None, 31, None),
        (20, 50, None, None, 1013.25),
        (20, 50, None, None, None),
        (20, None, 0.2, 31, 1013.25),
        (20, None, 0.2, 31, None),
        (20, None, 0.2, None, 1013.25),
        (20, None, 0.2, None, None),
        (20, None, None, 31, 1013.25),
        (20, None, None, 31, None),
        (20, None, None, None, 1013.25),
        (20, None, None, None, None),
        (None, 50, 0.2, 31, 1013.25),
        (None, 50, 0.2, 31, None),
        (None, 50, 0.2, None, 1013.25),
        (None, 50, 0.2, None, None),
        (None, 50, None, 31, 1013.25),
        (None, 50, None, 31, None),
        (None, 50, None, None, 1013.25),
        (None, 50, None, None, None),
        (None, None, 0.2, 31, 1013.25),
        (None, None, 0.2, 31, None),
        (None, None, 0.2, None, 1013.25),
        (None, None, 0.2, None, None),
        (None, None, None, 31, 1013.25),
        (None, None, None, 31, None),
        (None, None, None, None, 1013.25),
        (None, None, None, None, None),
    ),
)
def test_pet_static_values_is_none(ta, rh, v, tmrt, p):
    with pytest.raises(TypeError) as exc_info:
        pet_static(ta=ta, rh=rh, v=v, tmrt=tmrt, p=p)
    assert "can't be converted to double" in exc_info.value.args[0]


@pytest.mark.parametrize('f', [mrt, mrt_np])
@pytest.mark.parametrize(
    ('tg', 'va', 'ta', 'd', 'expected'),
    (
        # the first two cases are example from DIN EN ISO 7726, however
        # they perform very strict rounding, leading to slightly different results
        # here since there is no rounding done
        pytest.param(55, 0.3, 30, 0.15, 74.4, id='forced convection'),
        pytest.param(53.2, 0.3, 30, 0.1, 74.6, id='globe diameter 0.1'),
        pytest.param(55, 0.1, 30, 0.15, 55.0, id='natural convection'),
    ),
)
def test_tmrt_scalar_values(tg, va, ta, d, expected, f):
    assert pytest.approx(f(tg=tg, va=va, ta=ta, d=d, e=0.95), abs=1e-1) == expected


def test_tmrt_array_values_default_d_e():
    tg = np.array([55, 53.2, 55])
    va = np.array([0.3, 0.3, 0.1])
    ta = np.array([30, 30, 30])
    expected = np.array([74.4, 71.6, 55.0])
    assert pytest.approx(mrt_np(tg=tg, va=va, ta=ta), abs=1e-1) == expected


@pytest.mark.parametrize('f', [mrt, mrt_np])
def test_tmrt_array_values(f):
    tg = np.array([55, 53.2, 55])
    va = np.array([0.3, 0.3, 0.1])
    ta = np.array([30, 30, 30])
    d = np.array([0.15, 0.1, 0.15])
    e = np.array([0.95, 0.95, 0.95])
    expected = np.array([74.4, 74.6, 55.0])
    assert_array_almost_equal(
        f(tg=tg, va=va, ta=ta, d=d, e=e),
        expected,
        decimal=1,
    )


@pytest.mark.parametrize(
    ('ta', 'rh', 'expected'),
    (
        # this case is from Stull et al. 2011
        pytest.param(20, 50, 13.7),
    ),
)
def test_twb_scalar_values(ta, rh, expected):
    assert pytest.approx(twb(ta=ta, rh=rh), abs=1e-1) == expected
