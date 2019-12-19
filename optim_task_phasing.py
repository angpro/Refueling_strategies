import solver_fun_lambert as sfl
import solver_fun_phasing as sfp
import numpy as np
from astropy import units as u

R_EARTH = 6378
t_maneuvers = [1 * 3600 * 24 * 30, 3 * 3600 * 24 * 30]    # 1, 3 months
altitude = 35800 - R_EARTH      # km
FILE_NAME = "maneuvers_phasing.txt"
file_maneuvers = open(FILE_NAME, 'w')
file_maneuvers.write("Maneuvers on GEO orbit with altitude = 35800, for SC on Russia, number of SC = 5" + '\n')

longitude1 = 20
longitude2 = 160
sc_num = 5

distance_sc_dergee = (longitude2 - longitude1)/sc_num
sc_longitude = np.linspace(distance_sc_dergee/2+longitude1, longitude2-distance_sc_dergee/2, sc_num) / 180 * np.pi * u.rad
distance_sc = distance_sc_dergee / 180 * np.pi * u.rad

orbit = sfl.orbit_from_altitude(altitude=altitude)
period = orbit.period.value
file_maneuvers.write("Period of circ orbit = " + str(period) +'\n')

orbits = [orbit]
for i in range(1, sc_num+1):
    orbits.append(orbit.propagate_to_anomaly(sc_longitude[i-1]))

file_maneuvers.close()
# =========================================================================
# SSC is the first one

orbit_ssc_first = orbit.propagate_to_anomaly(sc_longitude[0]-distance_sc)
print("This is period rbit_ssc_first")
print(orbit_ssc_first.period)

# -------------------------------------------------------------------------
# phasing - one and multi-revolution

file_maneuvers = open(FILE_NAME, 'a')
file_maneuvers.write("-------------------------------------" + '\n')
file_maneuvers.write("Type of maneuver: multi-revolution phasing maneuver" + '\n')

num_dfi_in_fi_list = [1, 2, 3, 4]

for num_dfi_in_fi in num_dfi_in_fi_list:
    sfp.phasing_man(t_maneuvers, distance_sc_dergee, orbit_ssc_first,
                    file_maneuvers=file_maneuvers, num_dfi_in_fi=num_dfi_in_fi)

# num_dfi_in_fi = num_dfi_in_fi_list[1]
# t_maneuvers = [1 * 3600 * 24 * 30]
# sfp.phasing_man(t_maneuvers, distance_sc_dergee, orbit_ssc_first, type_of_ssc_location=1,
#                 file_maneuvers=file_maneuvers, num_dfi_in_fi=num_dfi_in_fi, plot=True)

file_maneuvers.close()
