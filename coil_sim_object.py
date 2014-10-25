# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 18:37:00 2014

@author: zchmlk
"""

import pylab as P
from scipy.constants import pi

from coil_set_object import CoilSet
from plasma_coil_object import PlasmaCoil

fdtype = P.float64


class CoilSim(object):

    def __init__(self, coil_set=CoilSet(), plasma_coil=PlasmaCoil()):

        self.coil_set = coil_set
        self.plasma_coil = plasma_coil

    def calc_psi(self):

        # combine coil floops with plasma floop
        r_floops_total_list = list(self.coil_set.r_floops)
        r_floops_total_list.extend([self.plasma_coil.r_floop])
        z_floops_total_list = list(self.coil_set.z_floops)
        z_floops_total_list.extend([self.plasma_coil.z_floop])

        r_floops_total = P.array(r_floops_total_list)
        z_floops_total = P.array(z_floops_total_list)

        self.r_floops_total = r_floops_total
        self.z_floops_total = z_floops_total
        self.n_floops = len(r_floops_total)

        # calc coils GFs for floops
        psi_coils_floops, b_r, b_z = \
            self.coil_set.calc_fields(r_range=self.r_floops_total,
                                      z_range=self.z_floops_total)
        self.psi_coils_floops = psi_coils_floops

        # calc plasma coil GFs for floops
        r_arr, z_arr, psi_plasma_coil_floops, b_r, b_z, r_out, z_out = \
            self.plasma_coil.calc_fields(r_range=self.r_floops_total,
                                         z_range=self.z_floops_total)
        self.psi_plasma_coil_floops = psi_plasma_coil_floops

        # combine all GFs for floops
        total_psi_floops_list = list(psi_coils_floops)
        total_psi_floops_list.append(psi_plasma_coil_floops)
        total_psi_floops = P.array(total_psi_floops_list)
        self.total_psi_floops = total_psi_floops

        # make mutual inductance matrix
        M = P.zeros([self.n_floops, self.n_floops], dtype=fdtype)

        # fill up M (NOT RIGHT, /100 AND OFF BY A BIT)
        M[:, 0:self.n_floops] = total_psi_floops

#        self.M = M
        # for some reason we need the transpose of M made this way
        self.M = M.T

        # solve for coil currents to make 1 Wb on all flux loops
        i_floops = P.linalg.solve(self.M, P.ones(self.n_floops, dtype=fdtype))

        # make GFs
        GFs = P.zeros([self.n_floops, 100, 100])

        # calc GFs for coils
        psi_coils, b_r, b_z = \
            self.coil_set.calc_fields(r_range=P.array([1e-6, 4.0],
                                                      dtype=fdtype),
                                      z_range=P.array([-1.5, 1.5],
                                                      dtype=fdtype),
                                      r_dim=100, z_dim=100)
        self.psi_coils = psi_coils

        # calc GF for plasma coil
        r_arr, z_arr, psi_plasma_coil, b_r, b_z, r_out, z_out = \
            self.plasma_coil.calc_fields(r_range=P.array([1e-6, 4.0],
                                                         dtype=fdtype),
                                         z_range=P.array([-1.5, 1.5],
                                                         dtype=fdtype),
                                         r_dim=100, z_dim=100)
        self.psi_plasma_coil = psi_plasma_coil

        # fill up GFs
        GFs[0:self.coil_set.n_coils, ...] = psi_coils
        GFs[-1, ...] = psi_plasma_coil

        self.GFs = GFs

        # normalize to I_p at 3.2 MA
        i_fact = 3.2e6 / i_floops[-1]
        i_floops *= i_fact

        # calculate flux on mesh
        psi = 0.0
        for ii, i_val in enumerate(i_floops):
            psi += i_val * GFs[ii, ...]

        self.i_floops = i_floops
        self.i_fact = i_fact
        self.psi = psi

        return psi

    def plot(self):

        self.calc_psi()

        r_range = P.array([1e-6, 4.0], dtype=fdtype)
        z_range = P.array([-1.5, 1.5], dtype=fdtype)
        r_dim = 100
        z_dim = 100

        # resistivity of Cu (Ohm-m)
        eta_Cu = 16.8e-9

        # coil resistances
        R_b_coils = eta_Cu * 2 * pi * self.coil_set.r_coils / \
            self.coil_set.r_widths / self.coil_set.z_widths

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

        # meshes for contour plots
        r_arr, z_arr = P.meshgrid(P.linspace(r_range[0], r_range[1], r_dim),
                                  P.linspace(z_range[0], z_range[1], z_dim))

        def pcoils(r_coils=self.coil_set.r_coils,
                   r_widths=self.coil_set.r_widths,
                   z_coils=self.coil_set.z_coils,
                   z_widths=self.coil_set.z_widths, lw=2):
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

        l_inj = 0.67 - 0.5 * w_i

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
        a_bs = self.coil_set.r_coils[-3] - 0.5 * \
            self.coil_set.r_widths[-3] - (r_u + l_inj)

        a_b = 0.7

        #r_m = 0.5 / r_bs * (a_b ** 2 - a_bs ** 2 + \
        #    (r_bs + r_0) ** 2 - r_0 ** 2)
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

        self.plasma_coil.plot_plasma()

        levels = P.arange(1.0, psi_max + delta_psi, delta_psi) * self.i_fact
        P.contour(r_arr, z_arr, self.psi, levels=levels, colors='w')

        levels = P.arange(psi_min, 1.0 + delta_psi, delta_psi) * self.i_fact
        P.matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        P.contour(r_arr, z_arr, self.psi, levels=levels, colors='gray')

        # need to plot black contours on midplane coil
        #rsub = P.where(r_arr[0, :] >= 5.0)[0]
        #levels = P.arange(1.0, psi_max + delta_psi, delta_psi) * i_fact
        #P.contour(r_arr[:, rsub], z_arr[:, rsub], psi[:, rsub],
        #          levels=levels, colors='gray')

        P.plot(self.coil_set.r_floops, self.coil_set.z_floops, 'ro')
        #P.plot(-r_floops, z_floops, 'ro')

        #no_text = True
        no_text = False
        # whether to annotate in MA or MW
        show_MA = True

        if not no_text:
            # annotate coil powers
            for ii in xrange(self.coil_set.n_coils):
                if self.i_floops[ii] >= 0:
                    signum = '+'
                else:
                    signum = '-'
                if show_MA:
                    tdata = 1e-6 * self.i_floops[ii]
                else:
                    tdata = 1e-6 * self.i_floops[ii] ** 2 * R_b_coils[ii]
        #        P.text(r_b_coils[ii], z_b_coils[ii], '%s%3.2f' %
        #               (signum, 1e-6 * i_floops[ii] ** 2 * R_b_coils[ii]),
        #               ha='center', va='center', backgroundcolor='w')
                P.text(self.coil_set.r_coils[ii], self.coil_set.z_coils[ii],
                       '%s%3.2f' % (signum, tdata),
                       ha='center', va='center', backgroundcolor='w')

            if self.i_floops[-1] >= 0:
                signum = '+'
            else:
                signum = '-'

            P.text(r_0, 0.0, '%s%3.2f MA' % (signum, 1e-6 * self.i_floops[-1]),
                   va='center', backgroundcolor='w')

        pcoils()
        P.axvline()

        coil_power = (self.i_floops[:self.coil_set.n_coils] ** 2 *
                      R_b_coils).sum()

        coil_current = abs(self.i_floops[:-1]).sum()

        if not no_text:
            P.title(r'Total Coil Power = %4.2f MW, Current %4.2f MA-turns' %
                    (1e-6 * coil_power, 1e-6 * coil_current))

        print r'Total Coil Power = %4.2f MW, Current %4.2f MA-turns' % \
              (1e-6 * coil_power, 1e-6 * coil_current)

        self.r_arr = r_arr
        self.z_arr = z_arr
        self.levels = levels
