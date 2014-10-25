# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 23:50:03 2014

@author: zchmlk
"""

import pylab as P

import filaments as F

from t3dinp import t3dinp

fdtype = P.float64
idtype = P.int16


class PlasmaCoil(object):

    def __init__(self, r_floop=0.5, z_floop=0.0,
                 i_p_coil_filename='hitpops.05.txt',
                 tris_filename='hitpops.05.t3d'):

        self.r_floop = r_floop
        self.z_floop = z_floop

        # read equilibrium file
        i_p_coils = P.loadtxt(i_p_coil_filename, delimiter=',', dtype=fdtype)
        self.i_p_coils = i_p_coils

        r_p_coils_full = i_p_coils[:, 0]
        z_p_coils_full = i_p_coils[:, 1]
        # ??? what is this scale factor, something to do with mu_0 ???
        beta = i_p_coils[:, 3] * 6.28e7
        i_p_coils_full = i_p_coils[:, 2]

        self.r_p_coils_full = r_p_coils_full
        self.z_p_coils_full = z_p_coils_full
        self.beta = beta
        self.i_p_coils_full = i_p_coils_full

        # choose subset where current is not zero

        sub = P.where(i_p_coils_full != 0.0)

        r_p_coils = r_p_coils_full[sub]
        z_p_coils = z_p_coils_full[sub]
        i_p_coils = i_p_coils_full[sub]
        n_p_coils = len(r_p_coils)

        self.r_p_coils = r_p_coils
        self.z_p_coils = z_p_coils
        self.i_p_coils = i_p_coils
        self.n_p_coils = n_p_coils

        r_p_widths = P.ones(n_p_coils, dtype=fdtype) * 0.05
        z_p_widths = 1.0 * r_p_widths
        n_r_p_filaments = P.ones(n_p_coils, dtype=idtype)
        n_z_p_filaments = 1 * n_r_p_filaments

        self.r_p_widths = r_p_widths
        self.z_p_widths = z_p_widths
        self.n_r_p_filaments = n_r_p_filaments
        self.n_z_p_filaments = n_z_p_filaments

        # read in triangle unstructured mesh information
        rzt, tris, pt = t3dinp(tris_filename)

        self.rzt = rzt
        self.tris = tris
        self.pt = pt

    def calc_fields(self, r_range=[0., 1], z_range=[-1, 1],
                    r_dim=20, z_dim=20, do_b_r=False, do_b_z=False):

        self.r_range = r_range
        self.z_range = z_range
        self.r_dim = r_dim
        self.z_dim = z_dim
        self.do_b_r = do_b_r
        self.do_b_z = do_b_z

        self.mesh = r_range, z_range, r_dim, z_dim, do_b_r, do_b_z

        # !!! not sure if 'mesh' will ever appear in this context !!!
        if 'mesh' not in dir(self):
            self.set_mesh()

        if len(self.r_range) == 2:
            r_arr, z_arr = P.meshgrid(P.linspace(self.r_range[0],
                                                 self.r_range[1],
                                                 self.r_dim),
                                      P.linspace(self.z_range[0],
                                                 self.z_range[1],
                                                 self.z_dim))
        else:
            r_arr, z_arr = self.r_range, self.z_range

        psi = 0.
        if self.do_b_r:
            b_r = 0.
        if self.do_b_z:
            b_z = 0.
        r_out = []
        z_out = []

        for ii, i_coil in enumerate(self.i_p_coils):
            I_0 = i_coil / self.n_r_p_filaments[ii] / self.n_z_p_filaments[ii]
            for kk in xrange(1, self.n_r_p_filaments[ii] + 1):
                r = self.r_p_coils[ii] + 0.5 * self.r_p_widths[ii] * \
                    (1. - (2. * kk - 1) / self.n_r_p_filaments[ii])
            for jj in xrange(1, self.n_z_p_filaments[ii] + 1):
                z = self.z_p_coils[ii] + 0.5 * self.z_p_widths[ii] * \
                    (1. - (2. * jj - 1) / self.n_z_p_filaments[ii])
                r_out.append(r)
                z_out.append(z)
                psi += F.poloidal_flux_filament(r_arr, z_arr, r, z, I_0)
                if self.do_b_r:
                    b_r += F.b_r_filament(r_arr, z_arr, r, z, I_0)
                if self.do_b_z:
                    b_z += F.b_z_filament(r_arr, z_arr, r, z, I_0)

        if not self.do_b_r:
            b_r = 0. * psi
        if not self.do_b_z:
            b_z = 0. * psi

        self.r_arr = r_arr
        self.z_arr = z_arr
        self.psi = psi
        self.b_r = b_r
        self.b_z = b_z
        self.r_out = r_out
        self.z_out = z_out

        return r_arr, z_arr, psi, b_r, b_z, r_out, z_out

    def plot_plasma(self):

        P.tricontourf(self.rzt[:, 0], self.rzt[:, 1],
                      self.tris, self.beta, 1001, zorder=0)
        cticks = P.linspace(0.0, 0.2, 5)
        P.colorbar(ticks=cticks, format='%.2f')
        P.jet()
