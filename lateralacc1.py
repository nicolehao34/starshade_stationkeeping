from astropy import units as u
from astropy import constants as const
import astropy
import numpy as np

class quantity:
    """
    This class stores variables such as Sun mass m1 and position vector
    of the sun r1.
    """
    m1 = const.M_sun
    m2 = const.M_earth
    # Assuming the Earth and Sun move on circular orbits around their
    # barycenter with a constant distance of 1AU between each other,
    # find the distance between the barycenter of the Earth Sun system and the Sun

    # Sun is at 0, Earth is 1 AU away
    coord1 = np.array([0,0,0])* u.au
    coord2 = np.array([1,0,0]) * u.au

    COM = ((m1*coord1+m2*coord2)/(m1+m2))

    # If the barycenter is at the origin, and the Sun is offset
    # in the negatuve x direction from the barycenter
    # r1, position vector of the Sun in COM frame
    r1 = coord1 - COM

    # r2, position vector of the Earth in COM frame
    r2 = coord2 - COM

# Sun mass
print(const.M_sun)
# Earth mass
print(const.M_earth)

# astronomical unit, average Sun-Earth distance
print(const.au)

# check what happens when dividing
# the astropy constamt by a float(like 2)
result = const.M_sun/2
print(result)

# try converting astronomical units to kilometers
q1 = 10. * u.au
q1.to(u.km)
print(q1)

# check quantity and unit
print(q1.value)
print(q1.unit)

# Write out the expression for the inertial acceleration of a satellite at r due
# to the inverse square law gravity from the Sun and Earth at r1 and r2.
# assuming r is a known vector of the satellite

r = np.array([0,0,0])
a1 = (- const.G*quantity.m1*(r-quantity.r1)/(np.linalg.norm(r-quantity.r1))**3).to(u.km/u.s**2)
a2 = (- const.G*quantity.m2*(r-quantity.r2)/(np.linalg.norm(r-quantity.r2))**3).to(u.km/u.s**2)

print("a1 is" + repr(a1))
print(repr(type(a1)))

print("a2 is" + repr(a2))
print(repr(type(a2)))

# the inertial acceleration of satellite is a1+a2
a_s = a1+a2


def find_acc(r):
    """
    calculate the inertial acceleraton from Sun and Earth given a position
    vector multiplied by some astropy units.
    Precondition: r is an astropy Quantity with units
    """
    # assert type(r) == astropy.Quantity
    # read documentation on astropy types
    m1 = const.M_sun
    m2 = const.M_earth

	#TODO: r1 and m1 are referenced here and are global variables.
	#You should put this into a class and have them as member or class variables
	#That way other code doesn't accidentally modify them unintentionally!

    a1 = (- const.G*quantity.m1*(r-quantity.r1)/(np.linalg.norm(r-quantity.r1)**3)).to(u.km/u.s**2)
    a2 = (- const.G*quantity.m2*(r-quantity.r2)/(np.linalg.norm(r-quantity.r2)**3)).to(u.km/u.s**2)
    a_s = (a1 + a2)

    return a_s

find_acc(np.array([0,0,0])*u.au)

# finding lateral differential acceleration
rt = np.array([1,0,0])*u.au
rs = np.array([-4,0,0])*u.au

a_t = find_acc(rt)
a_s = find_acc(rs)

print(np.linalg.norm(rs-rt))
axial_diff_acc = (np.dot((a_s-a_t),(rs-rt))*(rs-rt))/(np.linalg.norm(rs-rt))**2


def diff_lateral_acc(rt, rs):
    """
    This function takes the position of the telescope rt, a relative position
    vector for the starshade relaitve to the telescope rs.
    Returns differential lateral accelration.

    Precondition: rt, rs are both astropy quantities with units.
    """
    # find acceleration of telescope and starshade
    a_t = find_acc(rt)
    a_s = find_acc(rs)

    #TODO: store a_s-a_t and rs-rt so that they don't need to be reomputed multiple times
    r_diff = rs - rt
    a_diff = a_s - a_t

    # calculate axial differential acceleration
    axial_diff_acc = (np.dot((a_diff),(r_diff))*(r_diff))/(np.linalg.norm(r_diff))**2

    # calculate lateral differential acceleration
    lateral_diff_acc = (a_diff) - axial_diff_acc

    return lateral_diff_acc


x = diff_lateral_acc(np.array([1,0,0])*u.au, np.array([-4,0,0])*u.au)
