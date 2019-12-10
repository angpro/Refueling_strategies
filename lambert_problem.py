from astropy import time
# epoch = time.Time("2015-05-09 10:43")  # UTC by default
from poliastro.bodies import Earth
from poliastro.twobody.orbit import Orbit
from poliastro.maneuver import Maneuver
from astropy import units as u
import numpy as np

# test
# ss0 = Orbit.from_body_ephem(Earth, date_launch)
# ssf = Orbit.from_body_ephem(Mars, date_arrival)
#man_lambert = Maneuver.lambert(ss0, ssf)

R_EARTH = 6378   # km
M = 3.98600*10**5   # km**3/s**2


def linear_velocity(radius_vector, semi_major_axis):
    lin_vel = (2*M*(1/radius_vector-1/(2*semi_major_axis)))**(1/2)
    return lin_vel  #float('{:.20f}'.format(lin_vel))


def period_from_altitude(altitude):    # meters
    altitude_tot = altitude + R_EARTH
    period = 2*np.pi*altitude_tot**(3/2)/M**0.5
    return period  # s


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
print("orbit1.a")
print(orbit1.a)
print("orbit1.period")
print(orbit1.period)

period_1 = period_from_altitude(altitude_1)

# time_maneuver_2 = float('{:.3f}'.format(2/3 * period_1 + period_1))
# time_maneuver_2 = float('{:.3f}'.format(2/3 * period_1 + period_1))
# time_maneuver_2 = 2/3 * period_1
time_maneuver_2 = float('{:.3f}'.format(period_1))
print("time_maneuver")
print(time_maneuver_2)


# orbit 2
# date_arrival = time.Time('2011-11-06 05:17', scale='tdb')
print("date_launch.unix")
print(date_launch.unix)
date_arrival = time.Time((date_launch.unix + time_maneuver_2), format='unix')
print("date_arrival in unix")
print(date_arrival)
date_arrival = date_arrival.iso
print("date_arrival in iso")
print(date_arrival)
# date_arrival = date_arrival[0: 16]
date_arrival = time.Time(date_arrival, scale='utc')
print("date_arrival time.Time(date_arrival, scale='utc'")
print(date_arrival)

# date_arrival = (date_launch.unix + time.Time(2/3 * time_maneuver_2)).iso
altitude_2 = altitude_1 # 15000.0  # km

linear_velocity_y_2 = float('{:.3f}'.format(linear_velocity(R_EARTH+altitude_2, R_EARTH+altitude_2)))

# r2_y = R_EARTH+altitude_2
alfa = 30. * np.pi / 180. # angle between x and r1 (sat1)

# r2_y = -(R_EARTH+altitude_2) * np.cos(alfa)
# r2_x = -(R_EARTH+altitude_2) * np.sin(alfa)

# r2 = [r2_x, r2_y, 0] * u.km
r2 = [-(R_EARTH+altitude_2) * np.sin(alfa), (R_EARTH+altitude_2) * np.cos(alfa), 0] * u.km
v2 = [0.000, linear_velocity_y_2, 0.000] * u.km / u.s
orbit2 = Orbit.from_vectors(Earth, r2, v2, epoch=date_arrival)


man_lambert = Maneuver.lambert(orbit1, orbit2, short=False)

print("Total cost")
print(man_lambert.get_total_cost())

print("Total time")
print(man_lambert.get_total_time())

dv_a, dv_b = man_lambert.impulses

print(dv_a)
print(dv_b)

print('r of orbit2')
print(orbit2.r)

print('epoch')
print(orbit2.epoch)
