import numpy as np
from scipy.special import elliprg


def ellipsoid(a, b, c, spacing=(1.0, 1.0, 1.0), levelset=False):
    """Generate ellipsoid for given semi-axis lengths.

    The respective semi-axis lengths are given along three dimensions in
    Cartesian coordinates. Each dimension may use a different grid spacing.

    Parameters
    ----------
    a : float
        Length of semi-axis along x-axis.
    b : float
        Length of semi-axis along y-axis.
    c : float
        Length of semi-axis along z-axis.
    spacing : 3-tuple of floats
        Grid spacing in three spatial dimensions.
    levelset : bool
        If True, returns the level set for this ellipsoid (signed level
        set about zero, with positive denoting interior) as np.float64.
        False returns a binarized version of said level set.

    Returns
    -------
    ellipsoid : (M, N, P) array
        Ellipsoid centered in a correctly sized array for given `spacing`.
        Boolean dtype unless `levelset=True`, in which case a float array is
        returned with the level set above 0.0 representing the ellipsoid.

    """
    if (a <= 0) or (b <= 0) or (c <= 0):
        raise ValueError('Parameters a, b, and c must all be > 0')

    offset = np.r_[1, 1, 1] * np.r_[spacing]

    # Calculate limits, and ensure output volume is odd & symmetric
    low = np.ceil(-np.r_[a, b, c] - offset)
    high = np.floor(np.r_[a, b, c] + offset + 1)

    for dim in range(3):
        if (high[dim] - low[dim]) % 2 == 0:
            low[dim] -= 1
        num = np.arange(low[dim], high[dim], spacing[dim])
        if 0 not in num:
            low[dim] -= np.max(num[num < 0])

    # Generate (anisotropic) spatial grid
    x, y, z = np.mgrid[
        low[0] : high[0] : spacing[0],
        low[1] : high[1] : spacing[1],
        low[2] : high[2] : spacing[2],
    ]

    if not levelset:
        arr = ((x / float(a)) ** 2 + (y / float(b)) ** 2 + (z / float(c)) ** 2) <= 1
    else:
        arr = ((x / float(a)) ** 2 + (y / float(b)) ** 2 + (z / float(c)) ** 2) - 1

    return arr


def ellipsoid_stats(a, b, c):
    """Calculate analytical volume and surface area of an ellipsoid.

    The surface area of an ellipsoid is given by

    .. math:: S=4\\pi b c R_G\\!\\left(1, \\frac{a^2}{b^2}, \\frac{a^2}{c^2}\\right)

    where :math:`R_G` is Carlson's completely symmetric elliptic integral of
    the second kind [1]_. The latter is implemented as
    :py:func:`scipy.special.elliprg`.

    Parameters
    ----------
    a : float
        Length of semi-axis along x-axis.
    b : float
        Length of semi-axis along y-axis.
    c : float
        Length of semi-axis along z-axis.

    Returns
    -------
    vol : float
        Calculated volume of ellipsoid.
    surf : float
        Calculated surface area of ellipsoid.

    References
    ----------
    .. [1] Paul Masson (2020). Surface Area of an Ellipsoid.
           https://analyticphysics.com/Mathematical%20Methods/Surface%20Area%20of%20an%20Ellipsoid.htm

    """
    if (a <= 0) or (b <= 0) or (c <= 0):
        raise ValueError('Parameters a, b, and c must all be > 0')

    # Volume
    vol = 4 / 3.0 * np.pi * a * b * c

    # Surface area
    surf = 3 * vol * elliprg(1 / a**2, 1 / b**2, 1 / c**2)

    return vol, surf
