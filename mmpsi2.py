import pylab as P

import filaments as F


def mmpsi(r_coils=[0.5, 0.5], r_widths=[0.05, 0.05], n_r_filaments=[10, 10],
          z_coils=[-0.25, 0.25], z_widths=[0.15, 0.15], n_z_filaments=[10, 10],
          i_coils=[1., 1.], r_range=[0., 1], z_range=[-1, 1],
          r_dim=20, z_dim=20, do_b_r=False, do_b_z=False):
    """
    'multiple make psi'
    calculates magnetic flux, and optionally radial and axial magnetic
    flux density for multiple coils at desired locations

    coils are input as sequences (array, list, tuple, etc.) with the
    following all having the same dimensions:

    r_coils         radial location of centroid of coil
    z_coils         axial location of centroid of coil
    r_widths        radial extent of coil
    z_widths        axial extent of coil
    n_r_filaments   number of filaments to approximate coil in radial direction
    n_z_filaments   number of filaments to approximate coil in axial direction
    i_coils         total current (Ampere-turns) in each coil

    flux (and fields) are calculated on the domain as given:

    r_range         if length is is 2 elements, it is assumed that r_range
    z_range         and z_range are each the desired min, max points for
    r_dim           r and z, with dimensions r_dim and z_dim
    z_dim

    if length of r_range is other than two, it is assumed desired locations
    for the calculation are in r_range and z_range, which must have the
    same shape (and r_dim and z_dim are ignored)
    """

    # if r_range has 2 elements, it is assumed that r_range and z_range
    # are each the desired min, max points for r and z, with dimensions
    # r_dim and z_dim
    if len(r_range) == 2:
        r_arr, z_arr = P.meshgrid(P.linspace(r_range[0], r_range[1], r_dim),
                                  P.linspace(z_range[0], z_range[1], z_dim))
    else:
        # otherwise, the user is passing in the r and z locations for
        # the desired calucation points
        r_arr, z_arr = r_range, z_range

    # poloidal flux (will become the same dimensions as r_arr and z_arr)
    psi = 0.
    # define variables for radial and axial field
    if do_b_r:
        b_r = 0.
    if do_b_z:
        b_z = 0.

    # list for r and z locations for all filaments
    r_filaments = []
    z_filaments = []

    # for each coil
    for ii, i_coil in enumerate(i_coils):
        # current in each filament
        I_0 = i_coil / n_r_filaments[ii] / n_z_filaments[ii]
        # over all radial filaments
        for kk in xrange(1, n_r_filaments[ii] + 1):
            # r location of filament
            r = r_coils[ii] + 0.5 * r_widths[ii] * \
                (1. - (2. * kk - 1) / n_r_filaments[ii])
            # over all axial filaments
            for jj in xrange(1, n_z_filaments[ii] + 1):
                # z location of filament
                z = z_coils[ii] + 0.5 * z_widths[ii] * \
                    (1. - (2. * jj - 1) / n_z_filaments[ii])
                # add r and z location to filaments
                r_filaments.append(r)
                z_filaments.append(z)
                # calculate flux over entire r_arr and z_arr for
                # filament at (r, z) with current I_0
                psi += F.poloidal_flux_filament(r_arr, z_arr, r, z, I_0)
                # likewise for B_R
                if do_b_r:
                    b_r += F.b_r_filament(r_arr, z_arr, r, z, I_0)
                # and B_Z
                if do_b_z:
                    b_z += F.b_z_filament(r_arr, z_arr, r, z, I_0)

    # dummy arrays for B_R and B_Z if they were NOT calculated
    if not do_b_r:
        b_r = 0. * psi
    if not do_b_z:
        b_z = 0. * psi

    # convert lists to NumPy arrays
    r_filaments, z_filaments = P.array(r_filaments), P.array(z_filaments)

    # return calculation domain, flux, B_R, B_Z, and filament locations
    return r_arr, z_arr, psi, b_r, b_z, r_filaments, z_filaments
