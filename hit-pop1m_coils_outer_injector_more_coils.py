# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 08:48:12 2012

@author: nelson
"""

from scipy.constants import pi
import pylab as P
from mmpsi import mmpsi
from t3dinp import t3dinp

# dtype definitions for arrays (change here will affect all)
fdtype = P.float64
idtype = P.int16

# resistivity of Cu (Ohm-m)
eta_Cu = 16.8e-9

# minor radius
a = 1.0

# major radius
r_0 = 1.5

# radius of vacuum vessel (from r_0), plus blanket if present
blanket_width = 0.0
r_vv0 = a + blanket_width

# radius of flux conserver to Uberplate feature
a_u = 0.242

# width of Uberplate
w_u = 0.35

# width of injector
w_i = 0.28

# r location of Uberplate (injector feature is perpendicular
#  to Uberplate and tangent to the FC)
r_u = r_0 + P.sqrt((2 * a - w_u) * a_u - w_u ** 2 / 4 + a ** 2)

# Alternative:
#  Could also specific height of Uberplate, z_u, and calculate a_u as:
#z_u = 1.95
#a_u = (w_u**2 / 4 + z_u**2 - a**2) / (2 * a - w_u)

# angle of mating point of FC and Uberplate feature
theta_u = P.arctan2(0.5 * w_u + a_u, a + a_u)

# center of upper FC to Uberplate feature arc
r_u_1 = r_u
z_u_1 = w_u / 2 + a_u

# center of lower FC to Uberplate feature arc
r_u_2 = r_u
z_u_2 = -(w_u / 2 + a_u)

# big coil parameters
b_coil_width = 0.1

b_coils_list = []
r_b_coils_list = []
z_b_coils_list = []
b_coils_r_width_list = []
b_coils_z_width_list = []
r_floop_list = []
z_floop_list = []

array_radius = 1.15

array_r_width = 0.1
array_z_width = 0.1

delta_theta_degrees = 15.0
delta_theta = P.deg2rad(delta_theta_degrees)

n_coils_each_side = 5

theta_start_degrees = 270.0
theta_start = P.deg2rad(theta_start_degrees)
theta_end = theta_start + (n_coils_each_side - 1) * delta_theta

for theta in P.linspace(theta_start, theta_end, n_coils_each_side):
    r_coil = r_0 + array_radius * P.cos(theta)
    z_coil = array_radius * P.sin(theta)
    r_b_coils_list.append(r_coil)
    z_b_coils_list.append(z_coil)
    r_floop = r_0 + a * P.cos(theta)
    z_floop = a * P.sin(theta)
    r_floop_list.append(r_floop)
    z_floop_list.append(z_floop)

theta_start_degrees = 90.0 - (n_coils_each_side - 1) * delta_theta_degrees
theta_start = P.deg2rad(theta_start_degrees)
theta_end = P.deg2rad(90.0)

for theta in P.linspace(theta_start, theta_end, n_coils_each_side):
    r_coil = r_0 + array_radius * P.cos(theta)
    z_coil = array_radius * P.sin(theta)
    r_b_coils_list.append(r_coil)
    z_b_coils_list.append(z_coil)
    r_floop = r_0 + a * P.cos(theta)
    z_floop = a * P.sin(theta)
    r_floop_list.append(r_floop)
    z_floop_list.append(z_floop)

# coils on cryoports
r_b_coils_list.extend([1.16 + 0.5 * b_coil_width, 1.16 + 0.5 * b_coil_width,
                       0.55 + b_coil_width, 0.55 + b_coil_width])
z_b_coils_list.extend([-1.1, 1.1, -1.1, 1.1])

# keep index for double-wide coils
i_wide_coil = len(z_b_coils_list) - 2

# RKT -> Change coils on injector here
#r_inj_coils = [2.9] * 2
#r_inj_coils = [2.8, 3.0, 3.2] * 2
#r_inj_coils = [2.9, 3.2] * 2
r_inj_coils = [2.9, 3.2, 3.5] * 2
n_inj_coils = len(r_inj_coils)

z_inj_coil = 0.25
z_inj_coils = P.ones(n_inj_coils, dtype=fdtype) * z_inj_coil

for ii in xrange(int(n_inj_coils / 2)):
    z_inj_coils[ii] *= -1.0

r_b_coils_list.extend(r_inj_coils)
z_b_coils_list.extend(z_inj_coils)

# coil at end of injector
inj_end_coil = False
if inj_end_coil:
    r_b_coils_list.append(1.5 + 1.17 + 0.67 + b_coil_width)
    z_b_coils_list.append(0.0)

# make arrays out of lists
r_b_coils = P.array(r_b_coils_list, dtype=fdtype)

z_b_upper = r_vv0 + 0.5 * b_coil_width

z_b_inj_delta = 0.1

z_b_coils = P.array(z_b_coils_list, dtype=fdtype)

# total number of big coils
n_b_coils = len(r_b_coils)

r_b_widths = P.ones(n_b_coils, dtype=fdtype) * b_coil_width
z_b_widths = 1.0 * r_b_widths

n_r_b_elements = P.ones(n_b_coils, dtype=idtype) * 10
n_z_b_elements = 1 * n_r_b_elements

# double-wide inner coils
r_b_widths[i_wide_coil:i_wide_coil + 2] *= 2.0
n_r_b_elements[i_wide_coil:i_wide_coil + 2] *= 2
#n_r_b_elements[-3:-1] *= 2

# list entry of small coil positions
r_s_width = 0.2
z_s_width = r_s_width

# RKT -> Change "s" coils here
r_mod = 0.05 * 0
z_mod = 0.05 * 0
r_s_coils = P.array([r_u - r_mod, r_u - r_mod], dtype=fdtype)
z_s_coils = P.array([0.5 * w_u + a_u - z_mod, -(0.5 * w_u + a_u - z_mod)],
                    dtype=fdtype)

# total number of small coils
n_s_coils = len(r_s_coils)

r_s_widths = P.ones(n_s_coils, dtype=fdtype) * r_s_width
z_s_widths = 1.0 * r_s_widths
n_r_s_elements = P.ones(n_s_coils, dtype=idtype) * 5
n_z_s_elements = 1 * n_r_s_elements

# plasma coils
i_p_coils = P.loadtxt('hitpops.05.txt', delimiter=',', dtype=fdtype)

r_p_coils_full = i_p_coils[:, 0]
z_p_coils_full = i_p_coils[:, 1]
beta = i_p_coils[:, 3] * 6.28e7
i_p_coils_full = i_p_coils[:, 2]

# choose subset where current is not zero
sub = P.where(i_p_coils_full != 0.0)

r_p_coils = r_p_coils_full[sub]
z_p_coils = z_p_coils_full[sub]
i_p_coils = i_p_coils_full[sub]

n_p_coils = len(r_p_coils)

r_p_widths = P.ones(n_p_coils, dtype=fdtype) * 0.05
z_p_widths = 1.0 * r_p_widths
n_r_p_elements = P.ones(n_p_coils, dtype=idtype)
n_z_p_elements = 1 * n_r_p_elements

# coil resistances
R_b_coils = eta_Cu * 2 * pi * r_b_coils / r_b_widths / z_b_widths
R_s_coils = eta_Cu * 2 * pi * r_s_coils / r_s_widths / z_s_widths
# fraction of Cu in cross-section
Cu_fract = 0.9
R_b_coils /= Cu_fract
R_s_coils /= Cu_fract

# flux loops
# common factor for big coil flux loops
fl_b_fact = 1.0 / P.sqrt((r_b_coils - r_0) ** 2 + z_b_coils ** 2)

# plasma filaments act as one coil
n_coils = n_b_coils + n_s_coils + 1

r_floops = P.zeros(n_coils, dtype=fdtype)
z_floops = 0.0 * r_floops

## big coil flux loops
# all on FC on line from (r_0, 0) to center of coil
r_floops[0:n_b_coils] = r_0 + a * (r_b_coils - r_0) * fl_b_fact
z_floops[0:n_b_coils] = a * z_b_coils * fl_b_fact

l_inj = 0.67 - 0.5 * w_i

# put last flux loop at end of injector
if inj_end_coil:
    for ii in xrange(n_b_coils - n_inj_coils - 1, n_b_coils - 1):
        r_floops[ii] = r_b_coils[ii]
        z_floops[ii] = z_b_coils[ii] / abs(z_b_coils[ii]) * 0.5 * w_i
    r_floops[n_b_coils - 1] = r_u + l_inj + 0.5 * w_i
else:
    for ii in xrange(n_b_coils - n_inj_coils, n_b_coils):
        r_floops[ii] = r_b_coils[ii]
        z_floops[ii] = z_b_coils[ii] / abs(z_b_coils[ii]) * 0.5 * w_i

# small coil flux floops
#  first two are 45 degrees into the arc transition to injectors
r_floops[n_b_coils] = r_u_1 + a_u * P.cos(1.25 * pi)
z_floops[n_b_coils] = z_u_1 + a_u * P.sin(1.25 * pi)

r_floops[n_b_coils + 1] = r_u_2 + a_u * P.cos(0.75 * pi)
z_floops[n_b_coils + 1] = z_u_2 + a_u * P.sin(0.75 * pi)

# plasma current flux loop
r_floops[-1] = r_0 - a
z_floops[-1] = 0.0

n_floops = len(r_floops)

# mutual inductance matrix
M = P.zeros([n_floops, n_floops], dtype=fdtype)

for ii in xrange(n_b_coils):
    # set up column vector for each coil's current
    # get flux from this coil on all flux loops
    r_arr, z_arr, psi, b_r, b_z, r_out, z_out = \
        mmpsi(r_coils=[r_b_coils[ii]], r_widths=[r_b_widths[ii]],
              n_r_elements=[n_r_b_elements[ii]],
              z_coils=[z_b_coils[ii]], z_widths=[z_b_widths[ii]],
              n_z_elements=[n_z_b_elements[ii]],
              i_coils=[1.0],
              r_range=r_floops, z_range=z_floops)
    # fill in mutual inductance column
    M[:, ii] = psi

for ii in xrange(n_s_coils):
    # set up column vector for each coil's current
    # get flux from this coil on all flux loops
    r_arr, z_arr, psi, b_r, b_z, r_out, z_out = \
        mmpsi(r_coils=[r_s_coils[ii]], r_widths=[r_s_widths[ii]],
              n_r_elements=[n_r_s_elements[ii]],
              z_coils=[z_s_coils[ii]], z_widths=[z_s_widths[ii]],
              n_z_elements=[n_z_s_elements[ii]],
              i_coils=[1.0],
              r_range=r_floops, z_range=z_floops)
    # fill in mutual inductance column
    M[:, n_b_coils + ii] = psi

if 'GF_fl_p' not in dir():
    print "Calculating plasma GF for flux loops, standby ..."
    r_arr, z_arr, psi, b_r, b_z, r_out, z_out = \
        mmpsi(r_coils=r_p_coils, r_widths=r_p_widths,
              n_r_elements=n_r_p_elements,
              z_coils=z_p_coils, z_widths=z_p_widths,
              n_z_elements=n_z_p_elements,
              i_coils=i_p_coils,
              r_range=r_floops, z_range=z_floops)
    print "Done"
    # fill in last column of mutual inductance
    GF_fl_p = psi

M[:, -1] = GF_fl_p

# solve for coil currents to make 1 Wb on all flux loops
i_floops = P.linalg.solve(M, P.ones(n_floops, dtype=fdtype))

# !!! fix i_floops here for if any coils are in series !!!

r_range = P.array([1e-6, 4.0], dtype=fdtype)
z_range = P.array([-1.5, 1.5], dtype=fdtype)
r_dim, z_dim = 100, 100

# Green's function matrix for mesh
GFs = P.zeros([n_floops, r_dim, z_dim])

if 'GFs_b_coils' not in dir():
    print "Calculating b_coils GF for mesh, standby ..."
    GFs_b_coils = P.zeros([n_b_coils, r_dim, z_dim])
    for ii in xrange(n_b_coils):
        r_arr, z_arr, psi, b_r, b_z, r_out, z_out = \
            mmpsi(r_coils=[r_b_coils[ii]], r_widths=[r_b_widths[ii]],
                  n_r_elements=[n_r_b_elements[ii]],
                  z_coils=[z_b_coils[ii]], z_widths=[z_b_widths[ii]],
                  n_z_elements=[n_z_b_elements[ii]],
                  i_coils=[1.0],
                  r_range=r_range, z_range=z_range,
                  r_dim=r_dim, z_dim=z_dim)
        GFs_b_coils[ii, ...] = psi
    print "Done"

if 'GFs_s_coils' not in dir():
    print "Calculating s_coils GF for mesh, standby ..."
    GFs_s_coils = P.zeros([n_s_coils, r_dim, z_dim])
    for ii in xrange(n_s_coils):
        r_arr, z_arr, psi, b_r, b_z, r_out, z_out = \
            mmpsi(r_coils=[r_s_coils[ii]], r_widths=[r_s_widths[ii]],
                  n_r_elements=[n_r_s_elements[ii]],
                  z_coils=[z_s_coils[ii]], z_widths=[z_s_widths[ii]],
                  n_z_elements=[n_z_s_elements[ii]],
                  i_coils=[1.0],
                  r_range=r_range, z_range=z_range,
                  r_dim=r_dim, z_dim=z_dim)
        GFs_s_coils[ii, ...] = psi
    print "Done"

if 'GF_p_coil' not in dir():
    print "Calculating plasma GF for mesh, standby ..."
    r_arr, z_arr, psi, b_r, b_z, r_out, z_out = \
        mmpsi(r_coils=r_p_coils, r_widths=r_p_widths,
              n_r_elements=n_r_p_elements,
              z_coils=z_p_coils, z_widths=z_p_widths,
              n_z_elements=n_z_p_elements,
              i_coils=i_p_coils,
              r_range=r_range, z_range=z_range, r_dim=r_dim, z_dim=z_dim)
    print "Done"
    GF_p_coil = psi

GFs[0:n_b_coils, ...] = GFs_b_coils
GFs[n_b_coils:n_b_coils + n_s_coils, ...] = GFs_s_coils
GFs[-1, ...] = GF_p_coil

# normalize to I_p at 3.2 MA
i_fact = 3.2e6 / i_floops[-1]
i_floops *= i_fact

# calculate flux on mesh
psi = 0.0
for ii, i_val in enumerate(i_floops):
    psi += i_val * GFs[ii, ...]

# meshes for contour plots
r_arr, z_arr = P.meshgrid(P.linspace(r_range[0], r_range[1], r_dim),
                          P.linspace(z_range[0], z_range[1], z_dim))

# lists of all coils for pcoils function
r_coils = r_b_coils.tolist()
r_coils.extend(r_s_coils.tolist())
z_coils = z_b_coils.tolist()
z_coils.extend(z_s_coils.tolist())
r_widths = r_b_widths.tolist()
r_widths.extend(r_s_widths.tolist())
z_widths = z_b_widths.tolist()
z_widths.extend(z_s_widths.tolist())


def pcoils(r_coils=r_coils, r_widths=r_widths,
           z_coils=z_coils, z_widths=z_widths, lw=2):
    # plots rectangular coils
    r_sq = P.array([0.5, -0.5, -0.5, 0.5, 0.5])
    z_sq = P.array([0.5, 0.5, -0.5, -0.5, 0.5])
    for r_coil, r_width, z_coil, z_width in zip(r_coils, r_widths,
                                                z_coils, z_widths):
        r = r_coil + r_sq * r_width
        z = z_coil + z_sq * z_width
        P.plot(r, z, color='b', lw=lw)
        # optional to plot coil on left side
#        if plot_left_side:
#            P.plot(-r, z, color='b')

# FC parameters
# inner midplane gap half-angle
delta_deg = 2
delta = delta_deg * pi / 180.0

# number of points for arcs, etc.
npts = 101

# FC arc from midplane gap to outside injector feature
theta_1 = P.linspace(pi + delta, 2 * pi - theta_u, npts)
theta_2 = P.linspace(theta_u, pi - delta, npts)

# outline of the FC, in two parts
r_FC_1 = r_0 + a * P.cos(theta_1)
z_FC_1 = a * P.sin(theta_1)
r_FC_2 = r_0 + a * P.cos(theta_2)
z_FC_2 = a * P.sin(theta_2)

# injector feature angles
theta_inj_1 = P.linspace(pi + theta_u, 1.5 * pi, npts)
theta_inj_2 = P.linspace(0.5 * pi, pi - theta_u, npts)

# outline of the injector feature, in two parts
r_inj_1 = r_u_1 + a_u * P.cos(theta_inj_1)
z_inj_1 = z_u_1 + a_u * P.sin(theta_inj_1)
r_inj_2 = r_u_2 + a_u * P.cos(theta_inj_2)
z_inj_2 = z_u_2 + a_u * P.sin(theta_inj_2)

#smaller curve feature from Uberplate to injector
r_tran = z_u_1 - a_u - 0.5 * w_i

theta_tran_1 = P.linspace(pi, 1.5 * pi, npts)
r_tran_1 = r_u + r_tran + r_tran * P.cos(theta_tran_1)
z_tran_1 = 0.5 * w_i + r_tran + r_tran * P.sin(theta_tran_1)

# default linewidth
lw = 2

# flux conserver color
FC_col = 'r'

P.clf()
P.axes(aspect='equal')
# plot FC sections and injector transition
P.plot(r_FC_1, z_FC_1, r_FC_2, z_FC_2, color=FC_col, lw=lw)
#P.plot(-r_FC_1, z_FC_1, -r_FC_2, z_FC_2, color=FC_col, lw=lw)
#P.plot(r_inj_1, z_inj_1, -r_inj_1, z_inj_1, color=FC_col, lw=lw)
#P.plot(r_inj_2, z_inj_2, -r_inj_2, z_inj_2, color=FC_col, lw=lw)
P.plot(r_inj_1, z_inj_1, color=FC_col, lw=lw)
P.plot(r_inj_2, z_inj_2, color=FC_col, lw=lw)

# injector to Uberplate transition feature
P.plot(r_tran_1, z_tran_1, color=FC_col, lw=lw)
P.plot(r_tran_1, -z_tran_1, color=FC_col, lw=lw)

# sides of injector
#P.plot([r_u, r_u + l_inj],
#       [0.5 * w_u, 0.5 * w_u], color=FC_col, lw=lw)
#P.plot([r_u, r_u + l_inj],
#       [-0.5 * w_u, -0.5 * w_u], color=FC_col, lw=lw)
P.plot([r_u + r_tran, r_u + l_inj],
       [0.5 * w_i, 0.5 * w_i], color=FC_col, lw=lw)
P.plot([r_u + r_tran, r_u + l_inj],
       [-0.5 * w_i, -0.5 * w_i], color=FC_col, lw=lw)
#P.plot([-(r_0 + 0.5 * w_u), -(r_0 + 0.5 * w_u)],
#       [z_u, z_u + 0.4], color=FC_col, lw=lw)
#P.plot([-(r_0 - 0.5 * w_u), -(r_0 - 0.5 * w_u)],
#       [z_u, z_u + 0.4], color=FC_col, lw=lw)

# end of injector
th = P.linspace(-0.5 * pi, 0.5 * pi, npts)
#rs = r_u + l_inj + 0.5 * w_u * P.cos(th)
rs = r_u + l_inj + 0.5 * w_i * P.cos(th)
#zs = 0.5 * w_u * P.sin(th)
zs = 0.5 * w_i * P.sin(th)
P.plot(rs, zs, color=FC_col, lw=lw)
#P.plot(-rs, zs, color=FC_col, lw=lw)

delta2_deg = 0.0
delta2 = delta2_deg * pi / 180.0

r_bs = r_u + l_inj
a_bs = r_b_coils[-1] - 0.5 * r_b_widths[-1] - (r_u + l_inj)

a_b = 0.7

#r_m = 0.5 / r_bs * (a_b ** 2 - a_bs ** 2 + (r_bs + r_0) ** 2 - r_0 ** 2)
#z_m = P.sqrt(a_b ** 2 - (r_m - r_0) ** 2)
#
#theta_m = P.arctan2(z_m, r_m)

#thm = P.linspace(-theta_m, theta_m, npts)

# shield for midplane coil
#delta_m = 17.4 * pi / 180.0
delta_m = -5.0 * pi / 180.0
thm = P.linspace(-0.5 * pi + delta_m, 0.5 * pi - delta_m, npts)

#P.plot(r_bs + a_bs * P.cos(thm), a_bs * P.sin(thm), color='k', lw=lw)

# angle to meet shield at midplane coil
#delta_b = 6.8 * pi / 180.0
delta_b = 11.5 * pi / 180.0

theta2 = P.linspace(-0.5 * pi, -delta_b, npts)
r_vv = r_0 + r_vv0 * P.cos(theta2)
z_vv = r_vv0 * P.sin(theta2)

#P.plot(r_vv, z_vv, color='k', lw=lw)
theta2 = P.linspace(delta_b, 0.5 * pi, npts)
r_vv = r_0 + r_vv0 * P.cos(theta2)
z_vv = r_vv0 * P.sin(theta2)

#P.plot(r_vv, z_vv, color='k', lw=lw)
#P.plot(-r_vv, z_vv, color='k', lw=lw)

#P.plot([0, r_vv[-1]], [z_vv[-1], z_vv[-1]], color='k', lw=lw)
#P.plot([0, r_vv[-1]], [-z_vv[-1], -z_vv[-1]], color='k', lw=lw)

#tripcolor(r_pr, z_pr, p_pr)

#levels = linspace(psi_mod.min(), psi_mod.max(), 30)
#levels = P.arange(-15.0, 15.0 + 0.5, 0.5) * i_fact
#delta_psi = 10.0
#psi_min = -20.0
#psi_max = 100.0
#levels = arange(psi_min, psi_max + delta_psi, delta_psi)

#P.contour(r_arr, z_arr, psi, levels=levels)
#P.prism()
#P.hsv()
#P.contour(-r_arr, z_arr, psi, levels=levels)
##P.prism()
#P.hsv()

delta_psi = 0.5
psi_min = -15.0
psi_max = 30.0

levels = P.arange(1.0, psi_max + delta_psi, delta_psi) * i_fact
P.contour(r_arr, z_arr, psi, levels=levels, colors='w')

levels = P.arange(psi_min, 1.0 + delta_psi, delta_psi) * i_fact
P.matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
P.contour(r_arr, z_arr, psi, levels=levels, colors='gray')

# need to plot black contours on midplane coil
#rsub = P.where(r_arr[0, :] >= 5.0)[0]
#levels = P.arange(1.0, psi_max + delta_psi, delta_psi) * i_fact
#P.contour(r_arr[:, rsub], z_arr[:, rsub], psi[:, rsub],
#          levels=levels, colors='gray')

P.plot(r_floops, z_floops, 'ro')
#P.plot(-r_floops, z_floops, 'ro')

if 'tris' not in dir():
    rzt, tris, pt = t3dinp('hitpops.05.t3d')

#do_conf = False
do_conf = True

if do_conf:
    P.tricontourf(rzt[:, 0], rzt[:, 1], tris, beta, 1001, zorder=0)
    cticks = P.linspace(0.0, 0.2, 5)
    P.colorbar(ticks=cticks, format='%.2f')
    P.jet()

#no_text = True
no_text = False
# whether to annotate in MA or MW
show_MA = True

if not no_text:
    # annotate coil powers
    for ii in xrange(n_b_coils):
        if i_floops[ii] >= 0:
            signum = '+'
        else:
            signum = '-'
        if show_MA:
            tdata = 1e-6 * i_floops[ii]
        else:
            tdata = 1e-6 * i_floops[ii] ** 2 * R_b_coils[ii]
#        P.text(r_b_coils[ii], z_b_coils[ii], '%s%3.2f' %
#               (signum, 1e-6 * i_floops[ii] ** 2 * R_b_coils[ii]),
#               ha='center', va='center', backgroundcolor='w')
        P.text(r_b_coils[ii], z_b_coils[ii], '%s%3.2f' %
               (signum, tdata),
               ha='center', va='center', backgroundcolor='w')

    if i_floops[-1] >= 0:
        signum = '+'
    else:
        signum = '-'

    P.text(r_0, 0.0, '%s%3.2f MA' % (signum, 1e-6 * i_floops[-1]),
           va='center', backgroundcolor='w')

    #for ii in [0, 2, 4]:
#    for ii in [0, 1]:
    for ii in xrange(n_s_coils):
        if i_floops[n_b_coils + ii] >= 0:
            signum = '+'
        else:
            signum = '-'
        if show_MA:
            tdata = 1e-6 * i_floops[n_b_coils + ii]
        else:
            tdata = 1e-6 * i_floops[n_b_coils + ii] ** 2 * \
                R_b_coils[n_b_coils + ii]
#        P.text(r_s_coils[ii], z_s_coils[ii],
#               '%s%3.2f' % (signum,
#                            1e-6 * i_floops[n_b_coils + ii] ** 2 *
#                            R_s_coils[ii]),
#               va='center', ha='center', backgroundcolor='w')
        P.text(r_s_coils[ii], z_s_coils[ii],
               '%s%3.2f' % (signum, tdata),
               va='center', ha='center', backgroundcolor='w')
    #for ii in [1, 3, 5]:
    #for ii in [1]:
    #    if i_floops[ii] >= 0:
    #        signum = '+'
    #    else:
    #        signum = '-'
    #    P.text(r_s_coils[ii] - 1.5 * r_s_widths[ii], z_s_coils[ii],
    #           '%s%3.2f' % (signum,
    #                        1e-6 * i_floops[n_b_coils + ii] ** 2 *
    #                        R_s_coils[ii]),
    #           va='center', ha='right', backgroundcolor='w')

    #text(r_s_coils[-1], z_s_coils[-1] + 1.5 * z_s_widths[-1], '%3.2f' %
    #     (1e-6 * R_s_coils[-1] * i_floops[n_b_coils + n_s_coils - 1]**2),
    #     ha='center', va='center', backgroundcolor='w')

pcoils()
P.axvline()

coil_power = (i_floops[:n_b_coils] ** 2 * R_b_coils).sum() + \
             (i_floops[n_b_coils:-1] ** 2 * R_s_coils).sum()

coil_current = abs(i_floops[:-1]).sum()

if not no_text:
    P.title(r'Total Coil Power = %4.2f MW, Current %4.2f MA-turns' %
            (1e-6 * coil_power, 1e-6 * coil_current))

print r'Total Coil Power = %4.2f MW, Current %4.2f MA-turns' % \
      (1e-6 * coil_power, 1e-6 * coil_current)