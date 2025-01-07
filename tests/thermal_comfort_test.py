import math

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_array_almost_equal
from pandas.testing import assert_series_equal
from thermal_comfort import heat_index
from thermal_comfort import mrt
from thermal_comfort import mrt_np
from thermal_comfort import pet_static
from thermal_comfort import twb
from thermal_comfort import utci_approx


def load_utci_test_data():
    with open('testing/utci_references.txt') as f:
        # skip the headers
        lines = f.readlines()[35:]
        values = []
        for line in lines:
            values.append([float(i) for i in line.split('\t')])

    return values


def load_heat_index_test_data():
    with open('testing/heat_index_reference.csv') as f:
        # skip the headers
        lines = f.readlines()[35:]
        values = []
        for line in lines:
            values.append([float(i) for i in line.split(',')])

    return values


@pytest.mark.filterwarnings('ignore:encountered a value for')
@pytest.mark.parametrize(
    (
        'ta', 'd_tmrt', 'va', 'rh', 'pa', 'offset',
        'utci', 'utci_table', 'utci_polynomial',
    ),
    (*load_utci_test_data(),),
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
    assert math.isnan(utci_approx(ta=ta, tmrt=tmrt, va=va, rh=rh))


@pytest.mark.filterwarnings('ignore:encountered a value for')
def test_utci_approx_with_vectors():
    data = np.array(load_utci_test_data())
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
        load_utci_test_data(),
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


@pytest.mark.xfail(reason='reshaping is not implemented yet')
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


@pytest.mark.parametrize('shape', [(2, 1), (2, 1, 1)])
def test_utci_approx_shapes_incorrect(shape):
    ta = np.array([20.5, 30.5]).reshape(shape)
    tmrt = np.array([50.5, 70.5]).reshape(shape)
    va = np.array([1.5, 2.5]).reshape(shape)
    rh = np.array([50.5, 60.5]).reshape(shape)
    with pytest.raises(TypeError) as excinfo:
        utci_approx(ta=ta, tmrt=tmrt, va=va, rh=rh)

    assert excinfo.value.args[0] == (
        'Only arrays with one dimension are allowed. '
        'Please reshape your array accordingly'
    )


def test_utci_approx_array_sizes_differ():
    ta = np.array([20.5])
    tmrt = np.array([50.5, 70.5])
    va = np.array([1.5, 2.5])
    rh = np.array([50.5, 60.5])
    with pytest.raises(ValueError) as excinfo:
        utci_approx(ta=ta, tmrt=tmrt, va=va, rh=rh)

    assert excinfo.value.args[0] == 'All arrays must have the same length'


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
def test_pet_static_missing_value(ta, rh, v, tmrt, p):
    assert math.isnan(pet_static(ta=ta, rh=rh, v=v, tmrt=tmrt, p=p))


def test_pet_static_missing_value_mixed_array():
    result = pet_static(
        ta=np.array([20, float('nan')]),
        rh=np.array([50, float('nan')]),
        v=np.array([0.5, float('nan')]),
        tmrt=np.array([20, float('nan')]),
        p=np.array([1013.5, float('nan')]),
    )
    assert_array_almost_equal(result, [18, float('nan')], decimal=1)


def test_pet_static_missing_value_mixed_series():
    df = pd.DataFrame({
        'ta': pd.Series([20, float('nan')]),
        'rh': pd.Series([50, float('nan')]),
        'v': pd.Series([0.5, float('nan')]),
        'tmrt': pd.Series([20, float('nan')]),
        'p': pd.Series([1013.5, float('nan')]),
    })
    df['result'] = pet_static(
        ta=df['ta'], rh=df['rh'],
        v=df['v'], tmrt=df['tmrt'], p=df['p'],
    )
    assert_series_equal(
        left=df['result'],
        right=pd.Series([18, float('nan')]),
        atol=1,
        check_names=False,
    )


@pytest.mark.parametrize('shape', [(2, 1), (2, 1, 1)])
def test_pet_static_shapes_incorrect(shape):
    ta = np.array([20.5, 30.5]).reshape(shape)
    tmrt = np.array([50.5, 70.5]).reshape(shape)
    v = np.array([1.5, 2.5]).reshape(shape)
    rh = np.array([50.5, 60.5]).reshape(shape)
    p = np.array([1013.5, 1010.5]).reshape(shape)
    with pytest.raises(TypeError) as excinfo:
        pet_static(ta=ta, rh=rh, v=v, tmrt=tmrt, p=p)

    assert excinfo.value.args[0] == (
        'Only arrays with one dimension are allowed. '
        'Please reshape your array accordingly'
    )


def test_pet_static_array_sizes_differ():
    ta = np.array([20.5])
    tmrt = np.array([50.5, 70.5])
    v = np.array([1.5, 2.5])
    rh = np.array([50.5, 60.5])
    p = np.array([1013.5, 1010.5])
    with pytest.raises(ValueError) as excinfo:
        pet_static(ta=ta, rh=rh, v=v, tmrt=tmrt, p=p)

    assert excinfo.value.args[0] == 'All arrays must have the same length'


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


@pytest.mark.parametrize('f', [mrt, mrt_np])
def test_tmrt_array_values_default_d_e(f):
    tg = np.array([55, 53.2, 55])
    va = np.array([0.3, 0.3, 0.1])
    ta = np.array([30, 30, 30])
    expected = np.array([74.4, 71.6, 55.0])
    assert pytest.approx(f(tg=tg, va=va, ta=ta), abs=1e-1) == expected


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


@pytest.mark.parametrize('f', [mrt, mrt_np])
@pytest.mark.parametrize('d', (0, -0.01, -100.5))
def test_mrt_d_outside_of_allowed_range(f, d):
    with pytest.raises(ValueError) as excinfo:
        f(tg=70, va=3, ta=20, d=d)

    assert excinfo.value.args[0] == 'The globe diameter (d) must be positive'


@pytest.mark.parametrize('f', [mrt, mrt_np])
@pytest.mark.parametrize('e', (-0.01, 1.1))
def test_mrt_e_outside_of_allowed_range(f, e):
    with pytest.raises(ValueError) as excinfo:
        f(tg=70, va=3, ta=20, e=e)

    assert excinfo.value.args[0] == 'The emissivity (e) must be between 0 and 1'


@pytest.mark.parametrize('shape', ((2, 1), (2, 1, 1)))
def test_mrt_shapes_incorrect(shape):
    tg = np.array([50.5, 70.5]).reshape(shape)
    ta = np.array([20.5, 30.5]).reshape(shape)
    va = np.array([1.5, 2.5]).reshape(shape)
    with pytest.raises(TypeError) as excinfo:
        mrt(tg=tg, va=va, ta=ta)

    assert excinfo.value.args[0] == (
        'Only arrays with one dimension are allowed. '
        'Please reshape your array accordingly'
    )


def test_mrt_array_sizes_differ():
    tg = np.array([50.5])
    ta = np.array([20.5, 30.5])
    va = np.array([1.5, 2.5])
    with pytest.raises(ValueError) as excinfo:
        mrt(tg=tg, va=va, ta=ta)

    assert excinfo.value.args[0] == 'All arrays must have the same length'


@pytest.mark.parametrize(
    ('ta', 'rh', 'expected'),
    (
        # this case is from Stull et al. 2011
        pytest.param(20, 50, 13.7),
    ),
)
def test_twb_scalar_values(ta, rh, expected):
    assert pytest.approx(twb(ta=ta, rh=rh), abs=1e-1) == expected


def test_twb_array_values():
    ta = np.array([20, 25])
    rh = np.array([50, 55])
    expected = np.array([13.7, 18.8])
    assert_array_almost_equal(twb(ta=ta, rh=rh), expected, decimal=1)


@pytest.mark.parametrize('shape', ((2, 1), (2, 1, 1)))
def test_twb_shapes_incorrect(shape):
    ta = np.array([20, 25]).reshape(shape)
    rh = np.array([50, 55]).reshape(shape)
    with pytest.raises(TypeError) as excinfo:
        twb(ta=ta, rh=rh)

    assert excinfo.value.args[0] == (
        'Only arrays with one dimension are allowed. '
        'Please reshape your array accordingly'
    )


def test_twb_array_sizes_differ():
    ta = np.array([20])
    rh = np.array([50, 55])
    with pytest.raises(ValueError) as excinfo:
        twb(ta=ta, rh=rh)

    assert excinfo.value.args[0] == 'All arrays must have the same length'


def _c2f(celsius):
    """convert celsius to fahrenheit"""
    return celsius * 9 / 5 + 32


def _f2c(fahrenheit):
    """convert fahrenheit to celsius"""
    return (fahrenheit - 32) * 5 / 9


@pytest.mark.parametrize(
    ('ta', 'rh', 'expected'),
    (*load_heat_index_test_data(),),
)
def test_heat_index_scalar_values(ta, rh, expected):
    assert _c2f(heat_index(ta=_f2c(ta), rh=rh)) == pytest.approx(expected, abs=1)


@pytest.mark.filterwarnings('ignore:encountered a value for')
def test_heat_index_array():
    data = np.array(load_heat_index_test_data())
    ta = data[:, 0]
    rh = data[:, 1]
    expected = data[:, 2]

    result = _c2f(heat_index(ta=_f2c(ta), rh=rh))
    assert_array_almost_equal(result, expected, decimal=0)


@pytest.mark.parametrize(
    ('ta', 'rh'),
    (
        pytest.param(26, 60, id='ta out of valid range'),
        pytest.param(40, 35, id='rh out of valid range'),
        pytest.param(10, 35, id='both out of valid range'),
    ),
)
def test_heat_index_scalar_values_out_of_range(ta, rh):
    assert math.isnan(heat_index(ta=ta, rh=rh))


@pytest.mark.parametrize('shape', ((2, 1), (2, 1, 1)))
def test_heat_index_shapes_incorrect(shape):
    ta = np.array([20, 25]).reshape(shape)
    rh = np.array([50, 55]).reshape(shape)
    with pytest.raises(TypeError) as excinfo:
        heat_index(ta=ta, rh=rh)

    assert excinfo.value.args[0] == (
        'Only arrays with one dimension are allowed. '
        'Please reshape your array accordingly'
    )


def test_heat_index_array_sizes_differ():
    ta = np.array([20])
    rh = np.array([50, 55])
    with pytest.raises(ValueError) as excinfo:
        heat_index(ta=ta, rh=rh)

    assert excinfo.value.args[0] == 'All arrays must have the same length'
