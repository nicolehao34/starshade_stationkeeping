import numpy as np
from matplotlib import pyplot as plt
import lateralacc1 as acc
from astropy import units as u

def to_cartesian(rho,theta,phi):
    """
    Returns the input vector in cartesian coordinates.

    This function converts a vector in spherical coordinates to cartesian coordinates.
    rho here is the magnitude of the r_rel (from Jackson's paper)
    phi is the azimuthal angle. theta is the polar angle.

    Returns cartesian coordinates using the arguments rho, theta, and phi.
    """
	# Convert to raidans
    phi = phi * np.pi/180.
    theta = theta * np.pi/180.

    x = rho*np.sin(phi)*np.cos(theta)
    y = rho*np.sin(phi)*np.sin(theta)
    z = rho*np.cos(phi)

    coord = np.array([x,y,z])
    return coord


def diff_lateral_acc_values(R, n, theta_values, phi_values):
    """
    This function creates a nested list of differential lateral acceleration
    values that correspond to the meshgrid of thetavalues and phivalues.

    n is the number of steps (dimension of the nested loop)
    thetavalues is a list of theta values that will be converted to cartesian
    coordinates and then input into the function diff_lateral_acc, which takes
    2 arguments, rt and rs. (Set rho as R, the magnitude of rs - rt
    when converting to cartesian coordinates.)

    Returns a nested list of corresponding lateral differential acceleration
    values of the same dimension as [theta, phi].
    """
    # create empty nested list of dimension nxn
    lst = [[]for i in range(n)]

    # fill nested list in with values of differential lateral acceleration
    for n in range(len(lst)):
        for i in range(len(lst)):
            # convert to cartesian coordinates
            r_rel = to_cartesian(R,theta_values[i],phi_values[n])

            # calculate differential lateral acceleration values
			#TODO: using a point somewhere around L2 (also make this an argument to the function) *******
            r_t = np.array([1.01,.001,.001])* u.au
            r_s = r_t + r_rel
            delta_a_l = acc.diff_lateral_acc(r_t,r_s* u.au)
            delta_a_l = np.linalg.norm(delta_a_l.value)

            # append value to list
            lst[n].append(delta_a_l)

    return lst


# create values for theta and phi
thetavalues = np.arange(0, 360, 5)
phivalues = np.arange(-90,90,2.5)


# create 2D grid
[theta, phi] = np.meshgrid(thetavalues, phivalues)


# subplot
fig, ax = plt.subplots(1,1) # subplot, 1 row and 1 column


# calculate differential lateral acceleration
delta_a_l = np.array(diff_lateral_acc_values(.00001,72, thetavalues, phivalues))


# create contour plot 
plt.contourf(theta, phi, delta_a_l,10) # units
plt.colorbar()
ax.set_title('Differential Lateral Acceleration (km/s^2)')
ax.set_xlabel('θ (deg)')
ax.set_ylabel('ϕ (deg)')

plt.show()
