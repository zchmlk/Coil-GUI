from scipy import sqrt
# uses the complete elliptic integral of the first kind,
# scipy.special.ellipk = integral(1/sqrt(1-m*sin(t)**2), t=0..pi/2)
# and the complete elliptic integral of the second kind,
# scipy.special.ellipe = integral(sqrt(1-m*sin(t)**2), t=0..pi/2)
from scipy.special import ellipk as K
from scipy.special import ellipe as E

# permeability of free space
from scipy.constants import mu_0
#mu_0 = 4e-7 * pi

mu_0_over_2_pi = 2e-7


def poloidal_flux_filament(r, z, r_coil, z_coil, I_0):
    """return poloidal flux (Wb) at location(s) (r, z) for a filament centered
    at (r_coil, z_coil) with current I_0
    """
    # parameters for flux formula
    d2 = (z - z_coil) * (z - z_coil) + (r + r_coil) * (r + r_coil)
    d = sqrt(d2)
    m = 4 * r * r_coil / d2

    # poloidal flux for upper inner coil
    return mu_0 * I_0 * d * ((1 - 0.5 * m) * K(m) - E(m))


def b_r_filament(r, z, r_coil, z_coil, I_0):

    # parameters for field formula
    d2 = (z - z_coil) * (z - z_coil) + (r + r_coil) * (r + r_coil)
    d = sqrt(d2)
    m = 4 * r * r_coil / d2

    return mu_0_over_2_pi * I_0 / d * (z - z_coil) / r * \
           ((r * r + r_coil * r_coil + (z - z_coil) * (z - z_coil)) /
            ((r - r_coil) * (r - r_coil) + (z - z_coil) *
             (z - z_coil)) * E(m) - K(m))


def b_z_filament(r, z, r_coil, z_coil, I_0):

    # parameters for field formula
    d2 = (z - z_coil) * (z - z_coil) + (r + r_coil) * (r + r_coil)
    d = sqrt(d2)
    m = 4 * r * r_coil / d2

    return mu_0_over_2_pi * I_0 / d * \
           (K(m) - (r * r - r_coil * r_coil + (z - z_coil) * (z - z_coil)) /
            ((r - r_coil) * (r - r_coil) +
             (z - z_coil) * (z - z_coil)) * E(m))