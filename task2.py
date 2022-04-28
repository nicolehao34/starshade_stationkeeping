# import numpy as np
# from matplotlib import pyplot as plt
# import lateralacc1 as acc
# from astropy import units as u
#
# # create values for theta and phi
# n
#
# # create 2D grid
# [theta, phi] = np.meshgrid(thetavalues, phivalues)
# fig, ax = plt.subplots(1,1) # subplot, 1 row and 1 column
#
#
# def to_cartesian(rho,theta,phi):
#     """
#     Returns the input vector in cartesian coordinates.
#
#     This function converts a vector in spherical coordinates to cartesian coordinates.
#     rho here is the magnitude of the r_rel (from Jackson's paper)
#     phi is the azimuthal angle. theta is the polar angle.
#
#     Returns cartesian coordinates using the arguments rho, theta, and phi.
#     """
#
#     x = rho*np.sin(phi)*np.cos(theta)
#     y = rho*np.sin(phi)*np.sin(theta)
#     z = rho*np.cos(phi)
#
#     coord = np.array([x,y,z])
#     return coord
#
#
# def diff_lateral_acc_values(n,theta_values, phi_values):
#     """
#     This function creates a nested list of differential lateral acceleration
#     values that correspond to the meshgrid of thetavalues and phivalues.
#
#     n is the number of steps (dimension of the nested loop)
#     thetavalues is a list of theta values that will be converted to cartesian
#     coordinates and then input into the function diff_lateral_acc, which takes
#     2 arguments, rt and rs. (Set rho as R, the magnitude of rs - rt
#     when converting to cartesian coordinates.)
#
#     Returns a nested list of corresponding lateral differential acceleration
#     values of the same dimension as [theta, phi].
#     """
#     # create empty nested list of dimension nxn
#     lst = [[]for i in range(n)]
#     # fill nested list in with values of differential lateral acceleration
#     for n in range(len(lst)):
#         for i in range(len(lst)):
#             R = np.sqrt((np.cos(theta_values[i])*np.cos(phi_values[n]))**2+(np.cos(theta_values[i])\
#                 *np.sin(phi_values[n]))**2+(np.sin(theta_values[i]))**2)
#
#             # convert to cartesian coordinates
#             r_rel = to_cartesian(R,theta_values[i],phi_values[n])
#
#             # calculate differential lateral acceleration values
#             delta_a_l = acc.diff_lateral_acc(np.array([0,0,0])* u.au,r_rel* u.au)
#             delta_a_l = np.linalg.norm(delta_a_l.value)
#
#             # append value to list
#             lst[n].append(delta_a_l)
#
#     # print(lst)
#     return lst
#
# delta_a_l = np.array(diff_lateral_acc_values(10, thetavalues, phivalues))
#
# plt.contourf(thetavalues, phivalues, delta_a_l, 10)
#
# # plt.plot(a_s, a_t, marker='.', color='k', linestyle='none')
#
# # ax.set_title('Filled Contour Plot')
# # ax.set_xlabel('polar angle theta')
# # ax.set_ylabel('azimuthal angle phi')
#
# plt.show()


import numpy as np
from matplotlib import pyplot as plt
import lateralacc1 as acc
from astropy import units as u
import pdb



def to_cartesian(rho,theta,phi):
    """
    Returns the input vector in cartesian coordinates.

    This function converts a vector in spherical coordinates to cartesian coordinates.
    rho here is the magnitude of the r_rel (from Jackson's paper)
    phi is the azimuthal angle. theta is the polar angle.

    Returns cartesian coordinates using the arguments rho, theta, and phi.
    """

    x = rho*np.sin(phi)*np.cos(theta)
    y = rho*np.sin(phi)*np.sin(theta)
    z = rho*np.cos(phi)

    coord = np.array([x,y,z])
    return coord


def diff_lateral_acc_values(n,theta_values, phi_values):
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
            R = np.sqrt((np.cos(theta_values[i])*np.cos(phi_values[n]))**2+(np.cos(theta_values[i])\
                *np.sin(phi_values[n]))**2+(np.sin(theta_values[i]))**2)

            # convert to cartesian coordinates
            r_rel = to_cartesian(R,theta_values[i],phi_values[n])

            # calculate differential lateral acceleration values
            delta_a_l = acc.diff_lateral_acc(np.array([0,0,0])* u.au,r_rel* u.au)
            delta_a_l = np.linalg.norm(delta_a_l.value)

            # append value to list
            lst[n].append(delta_a_l)

    # print(lst)
    return lst

# create values for theta and phi
thetavalues = np.arange(0, 180, 10)
phivalues = np.arange(-90,90,10)

# create 2D grid
[theta, phi] = np.meshgrid(thetavalues, phivalues)
fig, ax = plt.subplots(1,1) # subplot, 1 row and 1 column

delta_a_l = np.array(diff_lateral_acc_values(18, thetavalues, phivalues))

# pdb.set_trace()
plt.contourf(theta, phi, delta_a_l,10)

# plt.plot(a_s, a_t, marker='.', color='k', linestyle='none')

# ax.set_title('Filled Contour Plot')
# ax.set_xlabel('polar angle theta')
# ax.set_ylabel('azimuthal angle phi')

plt.show()
# plt.savefig("plot")
# pdb.set_trace()
