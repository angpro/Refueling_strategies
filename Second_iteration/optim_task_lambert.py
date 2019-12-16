import solver_fun_lambert as sfl
import solver_fun_phasing as sfp
import numpy as np
from astropy import units as u

R_EARTH = 6378
t_maneuvers = [1 * 3600 * 24 * 30, 3 * 3600 * 24 * 30]  # , 6 * 3600 * 24 * 30]   # 1, 3, 6 months
altitude = 35800 - R_EARTH  #km
FILE_NAME = "maneuvers_lambert.txt"
file_maneuvers = open(FILE_NAME, 'w')
file_maneuvers.write("Maneuvers on GEO orbit with altitude = 35800, for SC on Russia, number of SC = 5" + '\n')

longitude1 = 20
longitude2 = 160
sc_num = 5

distance_sc_dergee = (longitude2 - longitude1)/sc_num
sc_longitude = np.linspace(distance_sc_dergee/2+longitude1, longitude2-distance_sc_dergee/2, sc_num) / 180 * np.pi * u.rad
distance_sc = distance_sc_dergee / 180 * np.pi * u.rad

orbit = sfl.orbit_from_altitude(altitude=altitude)
orbits = [orbit]
for i in range(1, sc_num+1):
    orbits.append(orbit.propagate_to_anomaly(sc_longitude[i-1]))


# =========================================================================
# SSC is the first one

orbit_ssc_first = orbit.propagate_to_anomaly(sc_longitude[0]-distance_sc)

# -------------------------------------------------------------------------
# lambert - one revolution
# ------------
# SSC is the first one, type_of_ssc_location=1
file_maneuvers = open(FILE_NAME, 'a')
file_maneuvers.write("-------------------------------------" + '\n')
file_maneuvers.write("Type of maneuver: general lambert maneuver" + '\n')

# for t_maneuver in t_maneuvers:
#     file_maneuvers.write("-- -- -- -- -- --" + '\n')
#     file_maneuvers.write("Total time for maneuver: " + str(t_maneuver) + '\n')
#
#     num_revol_list = [0, 1, 2, 3, 4]
#     for num_revol in num_revol_list:
#         file_maneuvers.write("Number of revolution: " + str(num_revol)+'\n')
#
#         sfl.maneuver_lambert(orbit_ssc_first, orbits[1], orbits[5], t_maneuver,
#                              file_maneuvers=file_maneuvers, distance_sc_dergee=distance_sc_dergee,
#                              type_of_ssc_location=1, num_revol=1, plot=False)

sfl.maneuver_lambert(orbit_ssc_first, orbits[1], orbits[5], t_maneuvers[1],
                     file_maneuvers=file_maneuvers, distance_sc_dergee=distance_sc_dergee,
                     type_of_ssc_location=1, num_revol=6, plot=True)

file_maneuvers.close()
