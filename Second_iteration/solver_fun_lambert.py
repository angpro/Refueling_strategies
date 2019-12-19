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
ALTITUDE = 1150.0  # km


def check_alt(orbit, file_maneuvers):
    if orbit.a * (1-orbit.ecc) - R_EARTH * u.km < 300 * u.km:
        file_maneuvers.write("! This orbit LOWER THAN 300 km" + '\n')


def linear_velocity(radius_vector_mod, semi_major_axis):  # km
    lin_vel = (2*M*(1/radius_vector_mod-1/(2*semi_major_axis)))**(1/2)
    return lin_vel  #float('{:.20f}'.format(lin_vel))  # km/s


def orbit_from_altitude(altitude: ALTITUDE):  # km
    r_x = R_EARTH + altitude
    r = [r_x, 0, 0] * u.km
    linear_velocity_y = float('{:.3f}'.format(linear_velocity(R_EARTH + altitude, R_EARTH + altitude)))
    v = [0.000, linear_velocity_y, 0.000] * u.km / u.s

    orbit = Orbit.from_vectors(Earth, r, v, epoch=DATA_LAUNCH)
    return orbit


def orbit_from_orbit(init_orbit, part_of_period, alfa_init):    # alfa in degrees
    period_init = init_orbit.period
    t_maneuver = part_of_period * period_init
    mean_motion_init = init_orbit.n
    anomaly_init = (alfa_init/180 * np.pi + mean_motion_init * t_maneuver) * u.rad

    orbit_dprt = init_orbit.propagate_to_anomaly(anomaly_init)

    r = orbit_dprt.r
    v = orbit_dprt.v

    date_arrival = time.Time((DATA_LAUNCH.unix + t_maneuver), format='unix')

    orbit = Orbit.from_vectors(Earth, r, v, epoch=date_arrival)
    return orbit


def maneuver_lambert(orbit_dpt, orbit_arr_sc, orbit_last, t_maneuver, file_maneuvers, distance_sc_dergee, flag,
                     type_of_ssc_location=1, plot=False, num_revol=None,
                     is_short=False, intermediate=True):

    # calculate time to return back
    alfa_back = (360 - distance_sc_dergee*5)/360 * np.pi * u.rad
    orbit_back = orbit_last.propagate_to_anomaly(alfa_back)
    r2_back = orbit_back.r
    v2_back = orbit_back.v
    date_arrival_back = time.Time((DATA_LAUNCH.unix + 2.1 * orbit_dpt.period.value*((360 - distance_sc_dergee*5)/360)), format='unix')
    orbit_arr_back = Orbit.from_vectors(Earth, r2_back, v2_back, epoch=date_arrival_back)

    man_lambert_back = Maneuver.lambert(orbit_last, orbit_arr_back, short=is_short, M=1)

    cost_maneuver_back = man_lambert_back.get_total_cost()
    time_maneuver_back = man_lambert_back.get_total_time()

    dv_tot = cost_maneuver_back
    T_man_tot = time_maneuver_back

    if flag:
        file_maneuvers.write("Last maneuver to return at init point (similar for all), it is a part of maneuvers, "
                             "which is added automatic to final answers " +'\n')
        file_maneuvers.write("dv back: " + str(cost_maneuver_back) + '\n')
        file_maneuvers.write("T maneuver back: " + str(time_maneuver_back) + '\n\n')

    t_maneuver = (t_maneuver-time_maneuver_back.value)/5

    # if type_of_ssc_location == 1:
    #     file_maneuvers.write("Where is SSC: SSC is the first one or last one" + '\n')
    # elif type_of_ssc_location == 3:
    #     file_maneuvers.write("Where is SSC: SSC is in the middle" + '\n')

    file_maneuvers.write("Number of revolution to go to one sc: " + str(num_revol) + '\n')
    file_maneuvers.write("Time of maneuver to go to one sc [sec]: " + str(t_maneuver) + '\n')

    alfa = orbit_dpt.n.value * t_maneuver * u.rad
    orbit = orbit_arr_sc.propagate_to_anomaly(alfa)

    r2 = orbit.r
    v2 = orbit.v
    date_arrival = time.Time((DATA_LAUNCH.unix + t_maneuver), format='unix')
    orbit_arr = Orbit.from_vectors(Earth, r2, v2, epoch=date_arrival)

    # try:

    if num_revol:

        man_lambert = Maneuver.lambert(orbit_dpt, orbit_arr, short=is_short, M=num_revol)

    else:
        man_lambert = Maneuver.lambert(orbit_dpt, orbit_arr, short=is_short)

    # except Warning as e:
    #     print("WARNING!")
    #     print(e)
    #     file_maneuvers.write("! BAD condition" + '\n')
    #     return

    orbit_trans, orbit_target = orbit_dpt.apply_maneuver(man_lambert, intermediate=intermediate)

    if plot:
        # try to plot it
        fig, ax = plt.subplots()  # (figsize=(4, 4))
        plotter = StaticOrbitPlotter(ax)
        plotter.plot(orbit_dpt, label="Initial orbit", color="red")
        plotter.plot(orbit_trans, label="Trans orbit", color="blue")
        plotter.plot(orbit_target, label="Final orbit", color="green")
        if num_revol:
            plt.savefig("figures_lamb/ManeuverLambert, t_man [month] = " + str(t_maneuver/(3600 * 24 * 30)) + "T, revol=" + str(num_revol) + ".png")
        else:
            plt.savefig("figures_lamb/ManeuverLambert, t_man [month] = " + str(t_maneuver/(3600 * 24 * 30)) + "T.png")
        plt.show()

    check_alt(orbit_trans, file_maneuvers)

    cost_maneuver = man_lambert.get_total_cost()
    time_maneuver = man_lambert.get_total_time()

    dv_tot = dv_tot + cost_maneuver*5
    T_man_tot = T_man_tot + time_maneuver*5

    file_maneuvers.write("dv: " + str(dv_tot) + '\n')
    file_maneuvers.write("T maneuver: " + str(T_man_tot) + '\n\n')
    file_maneuvers.write("- - - - - - - - - - " + '\n')

    print("dv: " + str(dv_tot) + '\n')
    print("T maneuver: " + str(T_man_tot) + '\n\n')
    return False

    # -------------------------------------------------------
    # no go to init point
    # t_maneuver = t_maneuver/ 5
    #
    # if type_of_ssc_location == 1:
    #     file_maneuvers.write("Where is SSC: SSC is the first one or last one" + '\n')
    # elif type_of_ssc_location == 3:
    #     file_maneuvers.write("Where is SSC: SSC is in the middle" + '\n')
    #
    # file_maneuvers.write("Time of maneuver [month]: " + str(t_maneuver / (3600 * 24 * 30)) + '\n')
    #
    # alfa = orbit_dpt.n.value * t_maneuver * u.rad
    # orbit = orbit_arr_sc.propagate_to_anomaly(alfa)
    # r2 = orbit.r
    # v2 = orbit.v
    # date_arrival = time.Time((DATA_LAUNCH.unix + t_maneuver), format='unix')
    # orbit_arr = Orbit.from_vectors(Earth, r2, v2, epoch=date_arrival)
    #
    # if num_revol:
    #     man_lambert = Maneuver.lambert(orbit_dpt, orbit_arr, short=is_short, M=num_revol)
    # else:
    #     man_lambert = Maneuver.lambert(orbit_dpt, orbit_arr, short=is_short)
    #
    # cost_maneuver = man_lambert.get_total_cost()
    # time_maneuver = man_lambert.get_total_time()
    #
    # dv_tot_witout_ret = cost_maneuver * 6
    # T_man_tot_witout_ret = time_maneuver * 6
    #
    # file_maneuvers.write("dv without going in initial point: " + str(dv_tot_witout_ret) + '\n')
    # file_maneuvers.write("T maneuver without going in initial point: " + str(T_man_tot_witout_ret) + '\n')

    # return dv_tot, T_man_tot
    # dv_tot_witout_ret, T_man_tot_witout_ret
