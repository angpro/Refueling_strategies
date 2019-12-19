from astropy import time
# epoch = time.Time("2015-05-09 10:43")  # UTC by default
from poliastro.bodies import Earth
from poliastro.twobody.orbit import Orbit
from poliastro.maneuver import Maneuver
from astropy import units as u
import matplotlib.pyplot as plt
from poliastro.plotting import StaticOrbitPlotter
import numpy as np

R_EARTH = 6378   # km
M = 3.98600*10**5   # km**3/s**2
DATA_LAUNCH = time.Time('2011-10-26 15:02', scale='utc')
ALTITUDE = 35800 - R_EARTH


def linear_velocity(radius_vector_mod, semi_major_axis):  # km
    lin_vel = (2*M*(1/radius_vector_mod-1/(2*semi_major_axis)))**(1/2)
    return lin_vel  #float('{:.20f}'.format(lin_vel))  # km/s


r_x = R_EARTH + ALTITUDE
r = [r_x, 0, 0] * u.km
linear_velocity_y = float('{:.3f}'.format(linear_velocity(R_EARTH + ALTITUDE, R_EARTH + ALTITUDE)))
v = [0.000, linear_velocity_y, 0.000] * u.km / u.s

orbit_dpt = Orbit.from_vectors(Earth, r, v, epoch=DATA_LAUNCH)
period = orbit_dpt.period.value
t_maneuver = period + 28/360*period

alfa = (28/180 * np.pi + t_maneuver * orbit_dpt.n.value) * u.rad
orbit2 = orbit_dpt.propagate_to_anomaly(alfa)
r_2 = orbit2.r
v_2 = orbit2.v


print("Period: " + str(period))
print("t_maneuver = " + str(t_maneuver))

# 72699  # sec

date_arrival = time.Time((DATA_LAUNCH.unix + t_maneuver), format='unix')
r_arr = [r_x, 10, 0] * u.km
# orbit_arr = Orbit.from_vectors(Earth, r_2, v_2, epoch=date_arrival)

orbit_arr = Orbit.from_vectors(Earth, r_arr, v, epoch=date_arrival)

# orbit3 = orbit2.propagate_to_anomaly(2*30/180 * np.pi * u.rad)

man_lambert = Maneuver.lambert(orbit_dpt, orbit_arr, short=False, M=2)
orbit_trans, orbit_target = orbit_dpt.apply_maneuver(man_lambert, intermediate=True)

# try to plot it
fig, ax = plt.subplots()  # (figsize=(4, 4))
plotter = StaticOrbitPlotter(ax)
plotter.plot(orbit_dpt, label="Initial orbit", color="red")
plotter.plot(orbit_trans, label="Trans orbit", color="blue")
plotter.plot(orbit_target, label="Final orbit", color="green")

print("Total cost")
print(man_lambert.get_total_cost())

print("Total time")
print(man_lambert.get_total_time())

# plt.savefig(
#     "figures_lamb/ManeuverLambert, t_man [month] = " + str(t_maneuver / (3600 * 24 * 30)) + "T, revol=" + str(
#         num_revol) + ".png")
# plt.savefig("figures_lamb/ManeuverLambert, t_man [month] = " + str(t_maneuver / (3600 * 24 * 30)) + "T.png")
plt.show()
