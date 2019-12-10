from astropy import time
# epoch = time.Time("2015-05-09 10:43")  # UTC by default
from poliastro.bodies import Earth
from poliastro.twobody.orbit import Orbit
from poliastro.maneuver import Maneuver
from astropy import units as u
import matplotlib.pyplot as plt
from poliastro.plotting import StaticOrbitPlotter

# u.Quantity
import numpy as np

R_EARTH = 6378   # km
M = 3.98600*10**5   # km**3/s**2


def linear_velocity(radius_vector, semi_major_axis):  # km
    lin_vel = (2*M*(1/radius_vector-1/(2*semi_major_axis)))**(1/2)
    return lin_vel  #float('{:.20f}'.format(lin_vel))  # km/s


def period_from_altitude(altitude): # km
    altitude_tot = altitude + R_EARTH
    period = 2*np.pi*altitude_tot**(3/2)/M**0.5
    return period  # s


def mean_montion_from_altitude(altitude):  # km
    a = altitude + R_EARTH
    return (M/a**3)**0.5


# =========================================================================
# orbit 1
date_launch = time.Time('2011-10-26 15:02', scale='utc')
print("date_launch")
print(date_launch)
altitude_1 = 1150.0  # km

linear_velocity_y_1 = float('{:.3f}'.format(linear_velocity(R_EARTH+altitude_1, R_EARTH+altitude_1)))

r1_x = R_EARTH+altitude_1
r1 = [r1_x, 0, 0] * u.km
v1 = [0.000, linear_velocity_y_1, 0.000] * u.km / u.s
orbit1 = Orbit.from_vectors(Earth, r1, v1, epoch=date_launch)

t1 = period_from_altitude(altitude_1)

# =========================================================================
# =========================================================================
# orbit 2 for SC1

# CHANGE ME
# ------------------------
t_maneuver = 0.8 * t1 + t1
print("t_maneuver")
print(t_maneuver)
# dv = [0.1, 0.1, 0] * u.km / u.s
# ------------------------

alfa1 = 120 / 180 * np.pi + mean_montion_from_altitude(altitude_1) * t_maneuver

print("mean_montion_from_altitude(altitude_1) * t_maneuver")
print(mean_montion_from_altitude(altitude_1) * t_maneuver)

anomaly = alfa1 * u.rad

my_orbit1 = orbit1.propagate_to_anomaly(anomaly)
# my_orbit1 = orbit1.sample(values=1, min_anomaly=min_anomaly) # GCRS object

r2 = my_orbit1.r
v2 = my_orbit1.v #+ dv

date_arrival = time.Time((date_launch.unix + t_maneuver), format='unix')
print("date_arrival in iso")
print(date_arrival.iso)

orbit2 = Orbit.from_vectors(Earth, r2, v2, epoch=date_arrival)

# =========================================================================
# Maneuver1

man_lambert1 = Maneuver.lambert(orbit1, orbit2, short=False)
orbit_trans, orbit_target = orbit1.apply_maneuver(man_lambert1, intermediate=True)

'''
dv_a, dv_b = man_lambert1.impulses

print("This is maneuver Lambert 1")
print(dv_a)
print(dv_b)

dv_a_mod = (dv_a[1][0]**2+dv_a[1][1]**2+dv_a[1][2]**2)**0.5
dv_b_mod = (dv_b[1][0]**2+dv_b[1][1]**2+dv_b[1][2]**2)**0.5

dv_sum1 = dv_a_mod + dv_b_mod
print("This is dv_sum1")
print(dv_sum1)
'''

print("This is man_lambert1")
print(man_lambert1)

print("This is alfa1")
print(alfa1)
print(alfa1 / np.pi * 180)

print("Total cost")
print(man_lambert1.get_total_cost())

print("Total time")
print(man_lambert1.get_total_time())

# =========================================================================
# try to plot it

fig, ax = plt.subplots() #(figsize=(4, 4))
plotter = StaticOrbitPlotter(ax)

plotter.plot(orbit1, label="Initial orbit", color="red")
plotter.plot(orbit_trans, label="Trans orbit", color="blue")
plotter.plot(orbit_target, label="Final orbit", color="green")

plt.savefig("fig.png")
plt.show()
