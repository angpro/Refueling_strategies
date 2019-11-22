import numpy as np

# R_EARTH = 6378000    # meters
# M = 3.98600*10**14   # m**3/s**2

R_EARTH = 6378    # km
M = 3.98600*10**5   # km**3/s**2

def period_from_altitude(altitude: int):    # meters
    altitude_tot = altitude + R_EARTH
    period = 2*np.pi*altitude_tot**(3/2)/M**0.5
    return period  # s


def semi_major_axis_from_period(period: int):
    semi_major_axis = (M*period**2/(2*np.pi)**2)**(1/3)
    return semi_major_axis


def linear_velocity(radius_vector, semi_major_axis):
    lin_vel = (2*M*(1/radius_vector)-1/(2*semi_major_axis))**(1/2)
    return lin_vel  #float('{:.20f}'.format(lin_vel))


# task 1, symmetric

n_1 = 3
# altitude_1 = 1150000  # m
altitude_1 = 1150  # km

period_1 = period_from_altitude(altitude_1)
print("This is period [sec]: " + str(period_1))

time_maneuver = 2/3 * period_1 + period_1

semi_major_axis_maneuver_1 = semi_major_axis_from_period(time_maneuver)
print("This is R_EARTH+altitude_1 [km]: " + str(R_EARTH+altitude_1))
print("This is semi_major_axis_maneuver_1 [km]: " + str(semi_major_axis_maneuver_1))

print("Linear velocity1 for a = R_EARTH+altitude_1 [km/s]: " + str(linear_velocity(R_EARTH+altitude_1, R_EARTH+altitude_1)))

print("Linear velocity2 for a = semi_major_axis_maneuver_1 [km/s]: " + str(linear_velocity(R_EARTH+altitude_1, semi_major_axis_maneuver_1)))

delta_v_1 = linear_velocity(R_EARTH+altitude_1, R_EARTH+altitude_1) - linear_velocity(R_EARTH+altitude_1, semi_major_axis_maneuver_1)

print("This is delta_v for maneuver1 [km/s]: " + str(delta_v_1))
