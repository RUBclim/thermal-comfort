import pyperf


array_setup = '''\
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.random.seed(0)
SHAPE = 100_000
ta = np.random.uniform(20, 50, size=SHAPE)
tg = np.random.uniform(10, 80, size=SHAPE)
mrt_value = np.random.uniform(3, 80, size=SHAPE)
v = np.random.uniform(0, 15, size=SHAPE)
rh = np.random.uniform(10, 100, size=SHAPE)
p = np.random.uniform(980, 1050, size=SHAPE)
d = np.random.uniform(0.05, 0.15, size=SHAPE)
e = np.random.uniform(0.5, 1, size=SHAPE)
'''


def main() -> int:
    runner = pyperf.Runner()
    # TMRT
    runner.timeit(
        name='tmrt (numpy) scalar',
        stmt='mean_radiant_temp_np(ta=20, tg=50, v=3)',
        setup='from thermal_comfort import mean_radiant_temp_np',
    )
    runner.timeit(
        name='tmrt (numpy) array',
        stmt='mean_radiant_temp_np(ta=ta, tg=tg, v=v)',
        setup=f'from thermal_comfort import mean_radiant_temp_np\n{array_setup}',
    )
    runner.timeit(
        name='tmrt scalar',
        stmt='mean_radiant_temp(ta=20, tg=50, v=3, d=0.15, e=0.95)',
        setup='from thermal_comfort import mean_radiant_temp',
    )

    runner.timeit(
        name='tmrt array',
        stmt='mean_radiant_temp(ta=ta, tg=tg, v=v, d=d, e=e)',
        setup=f'from thermal_comfort import mean_radiant_temp\n{array_setup}',
    )

    # TWB
    runner.timeit(
        name='twb scalar',
        stmt='wet_bulb_temp(ta=20, rh=50)',
        setup='from thermal_comfort import wet_bulb_temp',
    )

    runner.timeit(
        name='twb array',
        stmt='wet_bulb_temp(ta=ta, rh=rh)',
        setup=f'from thermal_comfort import wet_bulb_temp\n{array_setup}',
    )

    # Heat Index
    runner.timeit(
        name='heat index scalar',
        stmt='heat_index(ta=20, rh=50)',
        setup='from thermal_comfort import heat_index',
    )

    runner.timeit(
        name='heat index array',
        stmt='heat_index(ta=ta, rh=rh)',
        setup=f'from thermal_comfort import heat_index\n{array_setup}',
    )

    runner.timeit(
        name='heat index extended scalar',
        stmt='heat_index_extended(ta=20, rh=50)',
        setup='from thermal_comfort import heat_index_extended',
    )

    runner.timeit(
        name='heat index extended array',
        stmt='heat_index_extended(ta=ta, rh=rh)',
        setup=f'from thermal_comfort import heat_index_extended\n{array_setup}',
    )
    # dew point
    runner.timeit(
        name='dew point scalar',
        stmt='dew_point(ta=20, rh=50)',
        setup='from thermal_comfort import dew_point',
    )

    runner.timeit(
        name='dew point array',
        stmt='dew_point(ta=ta, rh=rh)',
        setup=f'from thermal_comfort import dew_point\n{array_setup}',
    )

    # absolute humidity
    runner.timeit(
        name='absolute humidity scalar',
        stmt='absolute_humidity(ta=20, rh=50)',
        setup='from thermal_comfort import absolute_humidity',
    )

    runner.timeit(
        name='absolute humidity array',
        stmt='absolute_humidity(ta=ta, rh=rh)',
        setup=f'from thermal_comfort import absolute_humidity\n{array_setup}',
    )
    # specific humidity
    runner.timeit(
        name='specific humidity scalar',
        stmt='specific_humidity(ta=20, rh=50)',
        setup='from thermal_comfort import specific_humidity',
    )

    runner.timeit(
        name='specific humidity array',
        stmt='specific_humidity(ta=ta, rh=rh)',
        setup=f'from thermal_comfort import specific_humidity\n{array_setup}',
    )

    # UTCI
    runner.timeit(
        name='utci scalar',
        stmt='utci_approx(ta=20, tmrt=50, v=3, rh=50)',
        setup='from thermal_comfort import utci_approx',
    )

    runner.timeit(
        name='utci array',
        stmt='utci_approx(ta=ta, tmrt=mrt_value, v=v, rh=rh)',
        setup=f'from thermal_comfort import utci_approx\n{array_setup}',
    )

    # PET
    runner.timeit(
        name='pet scalar',
        stmt='pet_static(ta=20, tmrt=50, v=3, rh=50, p=1013.25)',
        setup='from thermal_comfort import pet_static',
    )

    runner.timeit(
        name='pet array',
        stmt='pet_static(ta=ta, tmrt=mrt_value, v=v, rh=rh, p=p)',
        setup=f'from thermal_comfort import pet_static\n{array_setup}',
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
