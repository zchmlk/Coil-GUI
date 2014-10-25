# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:57:15 2014

@author: zchmlk
"""

import pylab as P

import filaments as F


class Coil(object):
    """
    magnetic field Coil object

    set __init__() method for inputs

    Methods:

    calc_fields()   calculates flux (and optionally field) on desired domain
    """

    def __init__(self, i_coil=1.0, r_width=0.1, z_width=0.1, r_coil=0.5,
                 z_coil=0.25, n_r_filaments=6, n_z_filaments=6,
                 r_floop=1, z_floop=1):
        """
        Inputs:

        i_coil          total current in coil [Ampere-turns]
        r_coil          radial position of coil [m]
        z_coil          axial position of coil [m]
        r_width         radial extent of coil [m]
        z_width         axial extent of coil [m]
        n_r_filaments   number of filaments to approximate coil radially
        n_z_filaments   number of filaments to approximate coil axially
        r_floop         radial flux loop
        z_floop         axial flux loop
        """
        self.i_coil = i_coil
        self.r_width = r_width
        self.z_width = z_width
        self.r_coil = r_coil
        self.z_coil = z_coil
        self.n_r_filaments = n_r_filaments
        self.n_z_filaments = n_z_filaments
        self.r_floop = r_floop
        self.z_floop = z_floop

    """reference mmpsi2.py for calc_fields comments"""
    def calc_fields(self, r_range=[0., 1.], z_range=[-1., 1.],
                    r_dim=20, z_dim=20, do_b_r=False, do_b_z=False):

        if len(r_range) == 2:
            r_arr, z_arr = P.meshgrid(P.linspace(r_range[0], r_range[1],
                                                 r_dim),
                                      P.linspace(z_range[0], z_range[1],
                                                 z_dim))
        else:
            r_arr, z_arr = r_range, z_range

        psi = 0.

        if do_b_r:
            b_r = 0.
        if do_b_z:
            b_z = 0.

        r_filaments = []
        z_filaments = []

        for kk in xrange(1, self.n_r_filaments + 1):
            r = self.r_coil + 0.5 * self.r_width * \
                (1. - (2. * kk - 1) / self.n_r_filaments)
            for jj in xrange(1, self.n_z_filaments + 1):
                z = self.z_coil + 0.5 * self.z_width * \
                    (1. - (2. * jj - 1) / self.n_z_filaments)
                r_filaments.append(r)
                z_filaments.append(z)
                psi += F.poloidal_flux_filament(r_arr, z_arr, r, z,
                                                self.i_coil /
                                                self.n_r_filaments /
                                                self.n_z_filaments)
                if do_b_r:
                    b_r += F.b_r_filament(r_arr, z_arr, r, z, self.i_coil)
                if do_b_z:
                    b_z += F.b_z_filament(r_arr, z_arr, r, z, self.i_coil)
        if not do_b_r:
            b_r = 0. * psi
        if not do_b_z:
            b_z = 0. * psi

        r_filaments, z_filaments = P.array(r_filaments), P.array(z_filaments)

        self.r_arr = r_arr
        self.z_arr = z_arr
        self.psi = psi
        self.b_r = b_r
        self.b_z = b_z
        self.r_filaments = r_filaments
        self.z_filaments = z_filaments

        return r_arr, z_arr, psi, b_r, b_z, r_filaments, z_filaments
