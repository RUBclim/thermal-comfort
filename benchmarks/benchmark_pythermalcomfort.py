import pyperf
from benchmark import array_setup


def main() -> int:
    runner = pyperf.Runner()
    # TMRT
    runner.timeit(
        name='tmrt scalar',
        stmt='t_mrt(tg=50, tdb=20, v=3, standard="ISO")',
        setup='from pythermalcomfort.psychrometrics import t_mrt',
    )

    runner.timeit(
        name='tmrt array',
        stmt='t_mrt(tg=tg, tdb=ta, v=va, standard="ISO")',
        setup=f'from pythermalcomfort.psychrometrics import t_mrt\n{array_setup}',
    )
    # TWB
    runner.timeit(
        name='twb scalar',
        stmt='t_wb(tdb=20, rh=50)',
        setup='from pythermalcomfort.psychrometrics import t_wb',
    )
    runner.timeit(
        name='twb array',
        stmt='t_wb(tdb=ta, rh=rh)',
        setup=f'from pythermalcomfort.psychrometrics import t_wb\n{array_setup}',
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
        stmt='f(tdb=ta, tr=mrt_value, v=va, rh=rh, p_atm=p, met=1.37, clo=0.5)',
        setup=(
            f'from pythermalcomfort.models import pet_steady\n{array_setup}'
            f'f = np.vectorize(pet_steady, otypes=[float], cache=True)'
        ),
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
