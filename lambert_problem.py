from astropy import time
# epoch = time.Time("2015-05-09 10:43")  # UTC by default
from poliastro.bodies import Earth
from poliastro.twobody.orbit import Orbit
from poliastro.maneuver import Maneuver
from astropy import units as u

# test
# ss0 = Orbit.from_body_ephem(Earth, date_launch)
# ssf = Orbit.from_body_ephem(Mars, date_arrival)
#man_lambert = Maneuver.lambert(ss0, ssf)

R_EARTH = 6378   # km
M = 3.98600*10**5   # km**3/s**2


def linear_velocity(radius_vector, semi_major_axis):
    lin_vel = (2*M*(1/radius_vector-1/(2*semi_major_axis)))**(1/2)
    return lin_vel  #float('{:.20f}'.format(lin_vel))


# orbit 1
date_launch = time.Time('2011-10-26 15:02', scale='tdb')
altitude_1 = 1150.0  # km
linear_velocity_y_1 = float('{:.3f}'.format(linear_velocity(R_EARTH+altitude_1, R_EARTH+altitude_1)))
r1_x = R_EARTH+altitude_1
r1 = [r1_x, 0, 0] * u.km
v1 = [0.000, linear_velocity_y_1, 0.000] * u.km / u.s
orbit1 = Orbit.from_vectors(Earth, r1, v1, epoch=date_launch)


# orbit 2
date_arrival = time.Time('2011-11-06 05:17', scale='tdb')
altitude_2 = 15000.0  # km
linear_velocity_y_2 = float('{:.3f}'.format(linear_velocity(R_EARTH+altitude_2, R_EARTH+altitude_2)))
r2_y = R_EARTH+altitude_2
r2 = [0, r2_y, 0] * u.km
v2 = [0.000, linear_velocity_y_2, 0.000] * u.km / u.s
orbit2 = Orbit.from_vectors(Earth, r2, v2, epoch=date_arrival)


man_lambert = Maneuver.lambert(orbit1, orbit2)

dv_a, dv_b = man_lambert.impulses

print(dv_a)
print(dv_b)