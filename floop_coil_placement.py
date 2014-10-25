# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 18:37:00 2014

@author: zchmlk
"""

from scipy.constants import pi
import pylab as P
import gspread

'''
Constructs coils and their corresponding flux loops
Writes each coil to the specified google spreadsheet with each
row representing a coil and each column representing a coil parameter
***gspread must be installed on your device in order to run this script***
'''

# dtype definitions for arrays (change here will affect all)
fdtype = P.float64
idtype = P.int16

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

# angle of mating point of FC and Uberplate feature
theta_u = P.arctan2(0.5 * w_u + a_u, a + a_u)

# center of upper FC to Uberplate feature arc
r_u_1 = r_u
z_u_1 = w_u / 2 + a_u

# center of lower FC to Uberplate feature arc
r_u_2 = r_u
z_u_2 = -(w_u / 2 + a_u)

# coil parameters
coil_width = 0.1

coils_list = []
r_coils_list = []
z_coils_list = []
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
    r_coils_list.append(r_coil)
    z_coils_list.append(z_coil)
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
    r_coils_list.append(r_coil)
    z_coils_list.append(z_coil)
    r_floop = r_0 + a * P.cos(theta)
    z_floop = a * P.sin(theta)
    r_floop_list.append(r_floop)
    z_floop_list.append(z_floop)

# coils on cryoports
r_coils_list.extend([1.16 + 0.5 * coil_width, 1.16 + 0.5 * coil_width,
                       0.55 + coil_width, 0.55 + coil_width])
z_coils_list.extend([-1.1, 1.1, -1.1, 1.1])

# keep index for double-wide coils
i_wide_coil = len(z_coils_list) - 2

# Change coils on injector here
r_inj_coils = [2.9, 3.2, 3.5] * 2
n_inj_coils = len(r_inj_coils)

z_inj_coil = 0.25
z_inj_coils = P.ones(n_inj_coils, dtype=fdtype) * z_inj_coil

for ii in xrange(int(n_inj_coils / 2)):
    z_inj_coils[ii] *= -1.0

r_coils_list.extend(r_inj_coils)
z_coils_list.extend(z_inj_coils)

# coil at end of injector
inj_end_coil = False
if inj_end_coil:
    r_coils_list.append(1.5 + 1.17 + 0.67 + coil_width)
    z_coils_list.append(0.0)

# list entry of arc coil positions
r_a_width = 0.2
z_a_width = r_a_width

# Change arc coil params here
r_mod = 0.05 * 0
z_mod = 0.05 * 0
r_a_coils = P.array([r_u - r_mod, r_u - r_mod], dtype=fdtype)
z_a_coils = P.array([0.5 * w_u + a_u - z_mod, -(0.5 * w_u + a_u - z_mod)],
                    dtype=fdtype)

# total number of arc coils
n_a_coils = len(r_a_coils)

# add arc coils to r_coils, z_coils
for ii in xrange(0, n_a_coils):
    r_coils_list.append(r_a_coils[ii])
    z_coils_list.append(z_a_coils[ii])

# make arrays out of lists
r_coils = P.array(r_coils_list, dtype=fdtype)

z_upper = r_vv0 + 0.5 * coil_width

z_inj_delta = 0.1

z_coils = P.array(z_coils_list, dtype=fdtype)

# total number of coils
n_coils = len(r_coils)

r_widths = P.ones(n_coils, dtype=fdtype) * coil_width
z_widths = 1.0 * r_widths

n_r_filaments = P.ones(n_coils, dtype=idtype) * 10
n_z_filaments = 1 * n_r_filaments

# double-wide inner coils
r_widths[i_wide_coil:i_wide_coil + 2] *= 2.0
n_r_filaments[i_wide_coil:i_wide_coil + 2] *= 2

# flux loops
# common factor for coil flux loops
fl_fact = 1.0 / P.sqrt((r_coils - r_0) ** 2 + z_coils ** 2)

# one floop for each coil
r_floops = P.zeros(n_coils, dtype=fdtype)
z_floops = 0.0 * r_floops

# coil flux loops
# all on FC on line from (r_0, 0) to center of coil
r_floops[0:n_coils] = r_0 + a * (r_coils - r_0) * fl_fact
z_floops[0:n_coils] = a * z_coils * fl_fact

l_inj = 0.67 - 0.5 * w_i

# put last flux loop at end of injector
if inj_end_coil:
    for ii in xrange(n_coils - n_inj_coils - 1, n_coils - 1):
        r_floops[ii] = r_coils[ii]
        z_floops[ii] = z_coils[ii] / abs(z_coils[ii]) * 0.5 * w_i
    r_floops[n_coils - 1] = r_u + l_inj + 0.5 * w_i
else:
    for ii in xrange(n_coils - n_inj_coils, n_coils):
        r_floops[ii] = r_coils[ii]
        z_floops[ii] = z_coils[ii] / abs(z_coils[ii]) * 0.5 * w_i
    
# arc coil flux floops
# first two are 45 degrees into the arc transition to injectors
r_floops[n_coils - 2] = r_u_1 + a_u * P.cos(1.25 * pi)
z_floops[n_coils - 2] = z_u_1 + a_u * P.sin(1.25 * pi)

r_floops[n_coils - 1] = r_u_2 + a_u * P.cos(0.75 * pi)
z_floops[n_coils - 1] = z_u_2 + a_u * P.sin(0.75 * pi)

n_floops = len(r_floops)

# writes coil params to google spreadsheet

#@author: NJ
gc = gspread.login('zchmlk@gmail.com', '75579632!')
spread_sheet = gc.open("coil_set_spreadsheet")
worksheet = spread_sheet.sheet1
for ii in xrange(2, worksheet.row_count + 1):    
    worksheet.update_cell(ii,2, r_coils[ii - 2]);
    worksheet.update_cell(ii,3, z_coils[ii - 2]);
    worksheet.update_cell(ii,4, r_widths[ii - 2]);
    worksheet.update_cell(ii,5, z_widths[ii - 2]);
    worksheet.update_cell(ii,6, n_r_filaments[ii - 2]);
    worksheet.update_cell(ii,7, n_z_filaments[ii - 2]);
    worksheet.update_cell(ii,8, r_floops[ii - 2]);
    worksheet.update_cell(ii,9, z_floops[ii - 2]);

  
