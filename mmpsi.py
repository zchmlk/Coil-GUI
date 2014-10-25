import pylab as P

import filaments as F


def mmpsi(r_coils=[0.5, 0.5], r_widths=[0.05, 0.05], n_r_elements=[10, 10],
          z_coils=[-0.25, 0.25], z_widths=[0.15, 0.15], n_z_elements=[10, 10],
          i_coils=[1., 1.], r_range=[0., 1], z_range=[-1, 1],
          r_dim=20, z_dim=20, do_b_r=False, do_b_z=False):

    if len(r_range) == 2:
        r_arr, z_arr = P.meshgrid(P.linspace(r_range[0], r_range[1], r_dim),
                                  P.linspace(z_range[0], z_range[1], z_dim))
    else:
        r_arr, z_arr = r_range, z_range

    psi = 0.
    if do_b_r:
        b_r = 0.
    if do_b_z:
        b_z = 0.
    r_out = []
    z_out = []

    for ii, i_coil in enumerate(i_coils):
        I_0 = i_coil / n_r_elements[ii] / n_z_elements[ii]
        for kk in xrange(1, n_r_elements[ii] + 1):
            r = r_coils[ii] + 0.5 * r_widths[ii] * \
                (1. - (2. * kk - 1) / n_r_elements[ii])
            for jj in xrange(1, n_z_elements[ii] + 1):
                z = z_coils[ii] + 0.5 * z_widths[ii] * \
                    (1. - (2. * jj - 1) / n_z_elements[ii])
                r_out.append(r)
                z_out.append(z)
                psi += F.poloidal_flux_filament(r_arr, z_arr, r, z, I_0)
                if do_b_r:
                    b_r += F.b_r_filament(r_arr, z_arr, r, z, I_0)
                if do_b_z:
                    b_z += F.b_z_filament(r_arr, z_arr, r, z, I_0)

    if not do_b_r:
        b_r = 0. * psi
    if not do_b_z:
        b_z = 0. * psi

    return r_arr, z_arr, psi, b_r, b_z, r_out, z_out