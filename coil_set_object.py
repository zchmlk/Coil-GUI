# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 14:24:19 2014

@author: nan
"""

from coil_object import Coil
#import gspread
import pylab as P

fdtype = P.float64


class CoilSet(object):
    """
    set of Coil() objects
.
    set __init__() method for inputs

    Methods:

    calc_fields()   calculates flux (and optionally field)
                    of set of coils on desired domain
    contourf()      plots filled contour of flux on domain
    """

    def __init__(self, coils=[]):
        """
        Inputs:

        coils         array of Coil() objects

        """

        self.coils = coils
        self.n_coils = len(coils)

    def add_Coil(self, coil, r_floops_list=[], z_floops_list=[],
                 r_coils_list=[], z_coils_list=[],
                 r_widths_list=[], z_widths_list=[]):
        if type(coil) != Coil:
            print "Error in add_Coil method, input not Coil object"
            return

        self.coils.append(coil)
        self.n_coils = len(self.coils)

        r_floops_list.append(coil.r_floop)
        z_floops_list.append(coil.z_floop)
        r_coils_list.append(coil.r_coil)
        z_coils_list.append(coil.z_coil)
        r_widths_list.append(coil.r_width)
        z_widths_list.append(coil.z_width)

#        r_floops, z_floops, r_coils, z_coils, r_widths, z_widths = \
#            P.array(r_floops_list), P.array(z_floops_list), \
#            P.array(r_coils_list), P.array(z_coils_list), \
#            P.array(r_widths_list), P.array(z_widths_list)

        self.r_floops = P.array(r_floops_list)
        self.z_floops = P.array(z_floops_list)
        self.r_coils = P.array(r_coils_list)
        self.z_coils = P.array(z_coils_list)
        self.r_widths = P.array(r_widths_list)
        self.z_widths = P.array(z_widths_list)

    def remove_Coil(self, coil):
        if type(coil) != Coil:
            print "Error in remove_Coil method, input not Coil object"
            return

        coilSetTemp = []
        self.n_coils = 0

        for coilSetCoil in self.coils:
            if (coilSetCoil.r_coil != coil.r_coil) and \
               (coilSetCoil.z_coil != coil.z_coil):
                coilSetTemp.append(coilSetCoil)
        self.coils = None
        self.coils = []

        for coilTemp in coilSetTemp:
            self.add_Coil(coilTemp)

    def read_cvs(self, filename='coil_set_spreadsheet',
                 username='hitsiexperiment@gmail.com', password='2vinjcity'):
        '''@author: NJ
        before using this method, install gspread in your device
        Do not include empty rows in google spread sheet'''

        r_floops_list = []
        z_floops_list = []
        r_coils_list = []
        z_coils_list = []
        r_widths_list = []
        z_widths_list = []

        print "Reading in coil parameters, standby..."

        # login to a Google account
        gc = gspread.login(username, password)

        # read in Google spreadsheet
        spread_sheet = gc.open(filename)
        worksheet = spread_sheet.sheet1
        values = worksheet.get_all_values()

        for ii in xrange(1, len(values)):
            coil = Coil()
            coil.i_coil = float(values[ii][0])
            coil.r_coil = float(values[ii][1])
            coil.z_coil = float(values[ii][2])
            coil.r_width = float(values[ii][3])
            coil.z_width = float(values[ii][4])
            coil.n_r_filaments = int(float(values[ii][5]))
            coil.n_z_filaments = int(float(values[ii][6]))
            coil.r_floop = float(float(values[ii][7]))
            coil.z_floop = float(float(values[ii][8]))
            self.coils.append(coil)
            r_floops_list.append(coil.r_floop)
            z_floops_list.append(coil.z_floop)
            r_coils_list.append(coil.r_coil)
            z_coils_list.append(coil.z_coil)
            r_widths_list.append(coil.r_width)
            z_widths_list.append(coil.z_width)

        self.n_coils = len(self.coils)

#        r_floops, z_floops, r_coils, z_coils, r_widths, z_widths = \
#            P.array(r_floops_list), P.array(z_floops_list), \
#            P.array(r_coils_list), P.array(z_coils_list), \
#            P.array(r_widths_list), P.array(z_widths_list)
#
#        self.r_floops = r_floops
#        self.z_floops = z_floops
#        self.r_coils = r_coils
#        self.z_coils = z_coils
#        self.r_widths = r_widths
#        self.z_widths = z_widths

        self.r_floops = P.array(r_floops_list)
        self.z_floops = P.array(z_floops_list)
        self.r_coils = P.array(r_coils_list)
        self.z_coils = P.array(z_coils_list)
        self.r_widths = P.array(r_widths_list)
        self.z_widths = P.array(z_widths_list)

    def calc_fields(self, r_range=[0., 1], z_range=[-1, 1],
                    r_dim=20, z_dim=20,
                    do_b_r=False, do_b_z=False):

        self.r_range = r_range
        self.z_range = z_range
        self.r_dim = r_dim
        self.z_dim = z_dim
        self.do_b_r = do_b_r
        self.do_b_z = do_b_z

        for coil in self.coils:
            coil.r_range = r_range
            coil.z_range = z_range
            coil.r_dim = r_dim
            coil.z_dim = z_dim
            coil.do_b_r = do_b_r
            coil.do_b_z = do_b_z

        self.mesh = r_range, z_range, r_dim, z_dim, do_b_r, do_b_z

        self.psi = 0.0
        self.b_r = 0.0
        self.b_z = 0.0

        psi_list = []
        b_r_list = []
        b_z_list = []

        for ii in xrange(self.n_coils):

            r_arr, z_arr, psi, b_r, b_z, r_filaments, z_filaments = \
                self.coils[ii].calc_fields(r_range=self.r_range,
                                           z_range=self.z_range,
                                           r_dim=self.r_dim,
                                           z_dim=self.z_dim,
                                           do_b_r=self.do_b_r,
                                           do_b_z=self.do_b_z)

            psi_list.append(psi)
            b_r_list.append(b_r)
            b_z_list.append(b_z)
            self.r_arr = r_arr
            self.z_arr = z_arr

        psi, b_r, b_z = P.array(psi_list), P.array(b_r_list), P.array(b_z_list)

        self.psi = psi
        self.b_r = b_r
        self.b_z = b_z

        psi_cumulative = 0.0

        for ii in xrange(self.n_coils):

            psi_cumulative += self.coils[ii].i_coil * self.psi[ii]

        self.psi_cumulative = psi_cumulative

        return psi, b_r, b_z
