from thermal_comfort import pet
from thermal_comfort import utci


def test_utci():
    res = pet(ta=20, RH=80, tmrt=70, v=3, mbody=1, age=2, ht=3, work=4, icl=5, sex=1)
    assert res == 0


def test_pet():
    res = utci(ta=20, rh=80, tmrt=70, va10m=3)
    assert res == 0
