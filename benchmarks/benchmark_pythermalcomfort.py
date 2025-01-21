import pyperf
from benchmark import array_setup


def main() -> int:
    runner = pyperf.Runner()
    # TMRT
    runner.timeit(
        name='tmrt scalar',
        stmt='mean_radiant_tmp(tg=50, tdb=20, v=3, standard="ISO")',
        setup='from pythermalcomfort.utilities import mean_radiant_tmp',
    )

    runner.timeit(
        name='tmrt array',
        stmt='mean_radiant_tmp(tg=tg, tdb=ta, v=va, standard="ISO")',
        setup=f'from pythermalcomfort.utilities import mean_radiant_tmp\n{array_setup}',
    )
    # TWB
    runner.timeit(
        name='twb scalar',
        stmt='wet_bulb_tmp(tdb=20, rh=50)',
        setup='from pythermalcomfort.utilities import wet_bulb_tmp',
    )
    runner.timeit(
        name='twb array',
        stmt='wet_bulb_tmp(tdb=ta, rh=rh)',
        setup=f'from pythermalcomfort.utilities import wet_bulb_tmp\n{array_setup}',
    )

    # Heat Index
    runner.timeit(
        name='heat index scalar',
        stmt='heat_index(tdb=20, rh=50)',
        setup='from pythermalcomfort.models import heat_index',
    )
    runner.timeit(
        name='heat index array',
        stmt='heat_index(tdb=ta, rh=rh)',
        setup=f'from pythermalcomfort.models import heat_index\n{array_setup}',
    )
    # dew point
    runner.timeit(
        name='dew point scalar',
        stmt='dew_point_tmp(tdb=20, rh=50)',
        setup='from pythermalcomfort.utilities import dew_point_tmp',
    )
    runner.timeit(
        name='dew point array',
        stmt='dew_point_tmp(tdb=ta, rh=rh)',
        setup=f'from pythermalcomfort.utilities import dew_point_tmp\n{array_setup}',
    )

    # UTCI
    runner.timeit(
        name='utci scalar',
        stmt='utci(tdb=20, tr=50, v=3, rh=50, limit_inputs=False)',
        setup='from pythermalcomfort.models import utci',
    )

    runner.timeit(
        name='utci array',
        stmt='utci(tdb=ta, tr=mrt_value, v=va, rh=rh, limit_inputs=False)',
        setup=f'from pythermalcomfort.models import utci\n{array_setup}',
    )

    # PET
    runner.timeit(
        name='pet scalar',
        stmt='pet_steady(tdb=20, tr=50, v=3, rh=50, p_atm=1013.25, met=1.37, clo=0.5)',
        setup='from pythermalcomfort.models import pet_steady',
    )

    runner.timeit(
        name='pet array',
        stmt='pet_steady(tdb=ta, tr=mrt_value, v=va, rh=rh, p_atm=p, met=1.37, clo=0.5)',  # noqa: E501
        setup=f'from pythermalcomfort.models import pet_steady\n{array_setup}',
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
