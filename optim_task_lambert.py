import solver_fun_lambert as sfl
import numpy as np
from astropy import units as u

R_EARTH = 6378
altitude = 35800 - R_EARTH  #km
FILE_NAME = "maneuvers_lambert.txt"
file_maneuvers = open(FILE_NAME, 'w')
file_maneuvers.write("Maneuvers on GEO orbit with altitude = 35800 km, for SC on Russia, number of SC = 5" + '\n')
file_maneuvers.close()

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
file_maneuvers.write("Type of maneuver: lambert maneuver" + '\n')
file_maneuvers.write("Where is SSC: SSC is the first one or last one" + '\n')

flag = True

t_man = 1.25 * 3600 * 24 * 30  # Change me
for num_revol in range(0, 7):
    file_maneuvers.write("Total and limit time of maneuver [sec]: " + str(t_man) + '\n')
    file_maneuvers.write("[month]: " + str(t_man/(3600 * 24 * 30)) + '\n')
    flag = sfl.maneuver_lambert(orbit_ssc_first, orbits[1], orbits[5], t_man,
                                file_maneuvers=file_maneuvers, distance_sc_dergee=distance_sc_dergee,
                                flag=flag, num_revol=num_revol, plot=True)

file_maneuvers.close()
