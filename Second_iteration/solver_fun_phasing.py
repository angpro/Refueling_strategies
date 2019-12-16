from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit
from astropy import units as u
from poliastro.bodies import Earth
import matplotlib.pyplot as plt
from poliastro.plotting import StaticOrbitPlotter
import numpy as np


ALTITUDE = 1150.0  # km
R_EARTH = 6378   # km
a = ALTITUDE + R_EARTH
M = 3.98600*10**5   # km**3/s**2
MAX_REVOL_ON_CIRC_ORBIT = 20


def check_alt(orbit):
    if orbit.a * (1-orbit.ecc) - R_EARTH * u.km < 300 * u.km:
        print("! This orbit lower then 300 km")
        file_maneuvers = open('maneuvers_phasing.txt', 'a')
        file_maneuvers.write("This orbit LOWER THAN 300 km" + '\n\n')
        print(orbit)
        return True


def period_from_altitude(altitude):    # km
    altitude_tot = altitude + R_EARTH
    period = 2*np.pi*altitude_tot**(3/2)/M**0.5
    return period  # s


def semi_major_axis_from_period(period):
    semi_major_axis = (M*period**2/(2*np.pi)**2)**(1/3)
    return semi_major_axis


def linear_velocity(radius_vector, semi_major_axis):
    lin_vel = (2*M*(1/radius_vector-1/(2*semi_major_axis)))**(1/2)
    return lin_vel  #float('{:.20f}'.format(lin_vel))


def maneuver_phasing_form_altitude(part_of_period, altitude=ALTITUDE, plot=False, intermediate=True):
    period = period_from_altitude(altitude)
    time_maneuver = part_of_period * period
    semi_major_axis_maneuver = semi_major_axis_from_period(time_maneuver)
    delta_v = linear_velocity(R_EARTH + altitude, R_EARTH + altitude) - \
              linear_velocity(R_EARTH + altitude, semi_major_axis_maneuver)

    if plot:
        a_ph = semi_major_axis_maneuver # + R_EARTH

        ss_i = Orbit.circular(Earth, alt=altitude * u.km)
        hoh = Maneuver.hohmann(ss_i, a_ph * u.km)
        ss_a, ss_f = ss_i.apply_maneuver(hoh, intermediate=intermediate)

        fig, ax = plt.subplots()  # (figsize=(4, 4))
        plotter = StaticOrbitPlotter(ax)
        plotter.plot(ss_i, label="Initial orbit", color="red")
        plotter.plot(ss_a, label="Phasing orbit", color="blue")
        plt.savefig("/figures/ManeuverPhasing, t_man=" + str(part_of_period) + ".png")

    return 2*delta_v, period*part_of_period


def maneuver_phasing_form_orbit(part_of_period, orbit, plot=False, intermediate=True):
    altitude = (orbit.a - R_EARTH * u.km).value
    period = (orbit.period).value
    time_maneuver = part_of_period * period
    semi_major_axis_maneuver = semi_major_axis_from_period(time_maneuver)
    delta_v = linear_velocity((orbit.a).value, (orbit.a).value) - \
              linear_velocity((orbit.a).value, semi_major_axis_maneuver)

    a_ph = semi_major_axis_maneuver # + R_EARTH

    ss_i = Orbit.circular(Earth, alt=altitude * u.km)
    hoh = Maneuver.hohmann(ss_i, a_ph * u.km)
    ss_a, ss_f = ss_i.apply_maneuver(hoh, intermediate=intermediate)

    if plot:
        fig, ax = plt.subplots()  # (figsize=(4, 4))
        plotter = StaticOrbitPlotter(ax)
        plotter.plot(ss_i, label="Initial orbit", color="red")
        plotter.plot(ss_a, label="Phasing orbit", color="blue")
        plt.savefig("figures/ManeuverPhasing, t_man=" + str(part_of_period) + ".png")

    return 2*delta_v, period*part_of_period, ss_a, ss_f


def phasing_man(t_maneuvers, distance_sc_dergee, orbit_ssc_first, type_of_ssc_location, file_maneuvers, num_dfi_in_fi=1):

    if type_of_ssc_location == 1:
        file_maneuvers.write("Where is SSC: SSC is the first one or last one" + '\n')
    elif type_of_ssc_location == 3:
        file_maneuvers.write("Where is SSC: SSC is in the middle" + '\n')
    file_maneuvers.write("Number of dfi in fi: " + str(num_dfi_in_fi) + '\n')

    for t_man in t_maneuvers:
        # t_man = t_man/5
        t_man = t_man/num_dfi_in_fi
        file_maneuvers.write("Time of maneuver [month] for one dfi: " + str(t_man / (3600 * 24 * 30)) + '\n')
        for num_revol_on_circ_orbit in range(0, MAX_REVOL_ON_CIRC_ORBIT):
            file_maneuvers.write("Number of revolutions on the circ orbit for one dfi: " + str(num_revol_on_circ_orbit) + '\n')

            if type_of_ssc_location == 1:
                dv, T_man, ss_a, ss_f = maneuver_phasing_form_orbit(
                part_of_period=num_revol_on_circ_orbit + distance_sc_dergee / 360,
                orbit=orbit_ssc_first, plot=False)

                if check_alt(ss_f):
                    return

                dv_tot_witout_ret = dv * 6 * num_dfi_in_fi
                T_man_tot_witout_ret = T_man * 6 * num_dfi_in_fi

                file_maneuvers.write("dv without going in initial point: " + str(dv_tot_witout_ret) + '\n')
                file_maneuvers.write("T maneuver without going in initial point: " + str(T_man_tot_witout_ret) + '\n')

                # return back
                dv, T_man, ss_a, ss_f = maneuver_phasing_form_orbit(
                    part_of_period=1 + (360 - 5 * distance_sc_dergee) / 360,
                    orbit=orbit_ssc_first, plot=False)

                dv_tot = dv_tot_witout_ret * 5 / 6 + dv
                T_man_tot = T_man_tot_witout_ret * 5 / 6 + T_man

                if t_man*num_dfi_in_fi < T_man_tot:
                    print("! OVER TIME")
                    file_maneuvers.write("OVER TIME" + '\n\n')
                    break

                file_maneuvers.write("dv: " + str(dv_tot) + '\n')
                file_maneuvers.write("T maneuver: " + str(T_man_tot) + '\n\n')

    maneuver_phasing_form_orbit(part_of_period=0 + distance_sc_dergee / 360,
                                                        orbit=orbit_ssc_first, plot=True)
    return dv_tot,T_man_tot, dv_tot_witout_ret, T_man_tot_witout_ret
