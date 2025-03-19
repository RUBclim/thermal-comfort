from thermal_comfort._thermal_comfort import thermal_comfort_mod
from thermal_comfort.models import heat_index
from thermal_comfort.models import heat_index_extended
from thermal_comfort.models import pet_static
from thermal_comfort.models import utci_approx
from thermal_comfort.utils import absolute_humidity
from thermal_comfort.utils import dew_point
from thermal_comfort.utils import mean_radiant_temp
from thermal_comfort.utils import sat_vap_press_ice
from thermal_comfort.utils import sat_vap_press_water
from thermal_comfort.utils import specific_humidity
from thermal_comfort.utils import wet_bulb_temp


# we expose the native fortran function as private functions if we want to ignore
# input validation. I.e. play stupid games, win stupid prizes, but get a speedup.
_utci_approx = thermal_comfort_mod.utci_approx
_pet_static = thermal_comfort_mod.pet_static
_mean_radiant_temp = thermal_comfort_mod.mean_radiant_temp
_wet_bulb_temp = thermal_comfort_mod.wet_bulb_temp
_heat_index = thermal_comfort_mod.heat_index
_heat_index_extended = thermal_comfort_mod.heat_index_extended
_sat_vap_press_water = thermal_comfort_mod.sat_vap_press_water
_sat_vap_press_ice = thermal_comfort_mod.sat_vap_press_ice
_dew_point = thermal_comfort_mod.dew_point
_absolute_humidity = thermal_comfort_mod.absolute_humidity
_specific_humidity = thermal_comfort_mod.specific_humidity

__all__ = [
    'utci_approx', '_utci_approx', 'pet_static',
    '_pet_static', 'mean_radiant_temp', '_mean_radiant_temp',
    'wet_bulb_temp', '_wet_bulb_temp', 'heat_index', '_heat_index',
    'heat_index_extended', '_heat_index_extended', 'sat_vap_press_water',
    '_sat_vap_press_water', 'sat_vap_press_ice', '_sat_vap_press_ice',
    'dew_point', '_dew_point', 'absolute_humidity', '_absolute_humidity',
    'specific_humidity', '_specific_humidity',
]
