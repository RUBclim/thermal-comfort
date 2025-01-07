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
va = np.random.uniform(0, 15, size=SHAPE)
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
        stmt='mrt_np(tg=50, va=3, ta=20)',
        setup='from thermal_comfort import mrt_np',
    )
    runner.timeit(
        name='tmrt (numpy) array',
        stmt='mrt_np(tg=tg, va=va, ta=ta)',
        setup=f'from thermal_comfort import mrt_np\n{array_setup}',
    )
    runner.timeit(
        name='tmrt scalar',
        stmt='mrt(tg=50, va=3, ta=20, d=0.15, e=0.95)',
        setup='from thermal_comfort import mrt',
    )

    runner.timeit(
        name='tmrt array',
        stmt='mrt(tg=tg, va=va, ta=ta, d=d, e=e)',
        setup=f'from thermal_comfort import mrt\n{array_setup}',
    )

    # TWB
    runner.timeit(
        name='twb scalar',
        stmt='twb(ta=20, rh=50)',
        setup='from thermal_comfort import twb',
    )

    runner.timeit(
        name='twb array',
        stmt='twb(ta=ta, rh=rh)',
        setup=f'from thermal_comfort import twb\n{array_setup}',
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

    # UTCI
    runner.timeit(
        name='utci scalar',
        stmt='utci_approx(ta=20, tmrt=50, va=3, rh=50)',
        setup='from thermal_comfort import utci_approx',
    )

    runner.timeit(
        name='utci array',
        stmt='utci_approx(ta=ta, tmrt=mrt_value, va=va, rh=rh)',
        setup=f'from thermal_comfort import utci_approx\n{array_setup}',
    )

    # PET
    runner.timeit(
        name='pet scalar',
        stmt='pet_static(ta=20, rh=50, v=3, tmrt=50, p=1013.25)',
        setup='from thermal_comfort import pet_static',
    )

    runner.timeit(
        name='pet array',
        stmt='pet_static(ta=ta, rh=rh, v=va, tmrt=mrt_value, p=p)',
        setup=f'from thermal_comfort import pet_static\n{array_setup}',
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
