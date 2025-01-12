from thermal_comfort._thermal_comfort import thermal_comfort_mod
from thermal_comfort.models import heat_index
from thermal_comfort.models import heat_index_extended
from thermal_comfort.models import pet_static
from thermal_comfort.models import utci_approx
from thermal_comfort.mrt import mrt
from thermal_comfort.mrt import mrt_np
from thermal_comfort.mrt import twb


# we expose the native fortran function as private functions if we want to ignore
# input validation. I.e. play stupid games, win stupid prizes, but get a speedup.
_utci_approx = thermal_comfort_mod.utci_approx
_pet_static = thermal_comfort_mod.pet_static
_mrt = thermal_comfort_mod.mrt
_twb = thermal_comfort_mod.twb
_heat_index = thermal_comfort_mod.heat_index
_heat_index_extended = thermal_comfort_mod.heat_index_extended

__all__ = [
    'utci_approx', '_utci_approx', 'pet_static',
    '_pet_static', 'mrt_np', 'mrt', '_mrt', 'twb', '_twb',
    'heat_index', '_heat_index', 'heat_index_extended', '_heat_index_extended',
]
