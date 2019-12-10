from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit
from astropy import units as u
from poliastro.bodies import Earth
import matplotlib.pyplot as plt
from poliastro.plotting import StaticOrbitPlotter


altitude_1 = 1150.0  # km
R_EARTH = 6378   # km
a = altitude_1 + R_EARTH

altitude_2_R = 9119.52 # - R_EARTH

# STEP1
ss_i_1 = Orbit.circular(Earth, alt=altitude_1 * u.km)
hoh_1 = Maneuver.hohmann(ss_i_1, altitude_2_R * u.km)

# ph = Maneuver.impulse([0.9836, 0, 0] * u.km / u.s)
#ss_a = ss_i.apply_maneuver(ph) #, intermediate=True)

print(hoh_1.get_total_cost())
print(hoh_1.get_total_cost())

ss_a_1, ss_f_1 = ss_i_1.apply_maneuver(hoh_1, intermediate=True)

# STEP2
hoh_2 = Maneuver.hohmann(ss_a_1, (altitude_1 + R_EARTH) * u.km)

ss_a_2, ss_f_2 = ss_i_1.apply_maneuver(hoh_1, intermediate=True)

'''
op = OrbitPlotter2D()
op.plot(orbit=ss_i, label="Initial orbit")
op.plot(orbit=ss_a, label="Transfer orbit")
op.plot(orbit=ss_f, label="Final orbit")
'''

fig, ax = plt.subplots() #(figsize=(4, 4))
plotter = StaticOrbitPlotter(ax)
# plotter = OrbitPlotter3D(figure=fig)
plotter.plot(ss_i_1, label="Initial orbit", color="red")
plotter.plot(ss_a_1, label="Phasing orbit", color="blue")
# plotter.plot(ss_f_2, label="Final orbit", color="green")

plt.savefig("fig.png")

plt.show()
