import pyperf
from benchmark import array_setup
reshape = '''\
mrt_value = mrt_value.reshape(-1, 2)
v = v.reshape(-1, 2)
rh = rh.reshape(-1, 2)
pa = 6.107 * (10. ** (7.5 * ta / (238. + ta)))
'''

person = '''\
from umep import PET_person
person = PET_person(mbody=75, age=35, height=1.75, activity=135, sex=1, clo=0.9)
'''


def main() -> int:
    runner = pyperf.Runner()

    # UTCI
    runner.timeit(
        name='utci scalar',
        stmt='utci_calculator(Ta=20, Tmrt=50, va10m=3, RH=50)',
        setup='from umep import utci_calculator',
    )

    runner.timeit(
        name='utci array',
        stmt='utci_calculator_grid(Ta=20, Tmrt=mrt_value, va10m=v, RH=50)',
        setup=f'from umep import utci_calculator_grid\n{array_setup}\n{reshape}',
    )

    # PET
    runner.timeit(
        name='pet scalar',
        stmt='_PET(ta=20, RH=50, tmrt=50, v=3, mbody=75, age=35, ht=1.75, work=135, icl=0.9, sex=1)',
        setup=f'from umep import _PET\n{person}',
    )

    runner.timeit(
        name='pet array',
        stmt='calculate_PET_grid(Ta=20, RH=50, Tmrt=mrt_value, va=v, pet=person)',  # noqa: E501
        setup=f'from umep import calculate_PET_grid\n{array_setup}\n{reshape}\n{person}',
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
