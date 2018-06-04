from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.stats
from scipy.special import i0


def von_mises_cdf_series(k,x,p):
    x = float(x)
    s = np.sin(x)
    c = np.cos(x)
    sn = np.sin(p*x)
    cn = np.cos(p*x)
    R = 0
    V = 0
    for n in range(p-1,0,-1):
        sn, cn = sn*c - cn*s, cn*c + sn*s
        R = 1./(2*n/k + R)
        V = R*(sn/n+V)

    return 0.5+x/(2*np.pi) + V/np.pi


def von_mises_cdf_normalapprox(k, x):
    b = np.sqrt(2/np.pi)*np.exp(k)/i0(k)
    z = b*np.sin(x/2.)
    return scipy.stats.norm.cdf(z)


def von_mises_cdf(k,x):
    ix = 2*np.pi*np.round(x/(2*np.pi))
    x = x-ix
    k = float(k)

    # These values should give 12 decimal digits
    CK = 50
    a = [28., 0.5, 100., 5.0]

    if k < CK:
        p = int(np.ceil(a[0]+a[1]*k-a[2]/(k+a[3])))

        F = np.clip(von_mises_cdf_series(k,x,p),0,1)
    else:
        F = von_mises_cdf_normalapprox(k, x)

    return F+ix
