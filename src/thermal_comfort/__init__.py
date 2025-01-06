from thermal_comfort._thermal_comfort import thermal_comfort_mod
from thermal_comfort.models import pet_static
from thermal_comfort.models import utci_approx
from thermal_comfort.mrt import mrt
from thermal_comfort.mrt import mrt_np


# we expose the native fortran function as private functions if we want to ignore
# input validation. I.e. play stupid games, win stupid prizes, but get a speedup.
_utci_approx = thermal_comfort_mod.utci_approx
_pet_static = thermal_comfort_mod.pet_static
_mrt = thermal_comfort_mod.mrt
twb = thermal_comfort_mod.twb

__all__ = [
    'utci_approx', '_utci_approx', 'pet_static',
    '_pet_static', 'mrt_np', 'mrt', '_mrt', 'twb',
]
