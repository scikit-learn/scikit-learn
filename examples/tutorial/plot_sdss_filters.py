"""
SDSS Filters
------------

This example downloads and plots the filters from the Sloan Digital Sky
Survey, along with a reference spectrum.
"""
import os
import urllib2

import numpy as np
import pylab as pl
from matplotlib.patches import Arrow

from sklearn.datasets import get_data_home

REFSPEC_URL = 'ftp://ftp.stsci.edu/cdbs/current_calspec/1732526_nic_002.ascii'
URL = 'http://www.sdss.org/dr7/instruments/imager/filters/%s.dat'

def fetch_filter(filter, data_home=None):
    data_home = get_data_home(data_home)
    assert filter in 'ugriz'
    url = URL % filter
    loc = os.path.join(data_home, '%s.dat' % filter)
    if not os.path.exists(loc):
        print "downloading from %s" % url
        F = urllib2.urlopen(url)
        open(loc, 'w').write(F.read())

    F = open(loc)
        
    data = np.loadtxt(F)
    return data

def fetch_vega_spectrum(data_home=None):
    data_home = get_data_home(data_home)
    refspec_file = os.path.join(data_home, REFSPEC_URL.split('/')[-1])
    if  not os.path.exists(refspec_file):
        print "downnloading from %s" % REFSPEC_URL
        F = urllib2.urlopen(REFSPEC_URL)
        open(refspec_file, 'w').write(F.read())

    F = open(refspec_file)

    data = np.loadtxt(F)
    return data


Xref = fetch_vega_spectrum()
Xref[:, 1] /= 2.1 * Xref[:, 1].max()

#----------------------------------------------------------------------
# Plot filters in color with a single spectrum
pl.figure()
pl.plot(Xref[:, 0], Xref[:, 1], '-k', lw=2)

for f,c in zip('ugriz', 'bgrmk'):
    X = fetch_filter(f)
    pl.fill(X[:, 0], X[:, 1], ec=c, fc=c, alpha=0.4)

kwargs = dict(fontsize=20, ha='center', va='center', alpha=0.5)
pl.text(3500, 0.02, 'u', color='b', **kwargs)
pl.text(4600, 0.02, 'g', color='g', **kwargs)
pl.text(6100, 0.02, 'r', color='r', **kwargs)
pl.text(7500, 0.02, 'i', color='m', **kwargs)
pl.text(8800, 0.02, 'z', color='k', **kwargs)

pl.xlim(3000, 11000)

pl.title('SDSS Filters and Reference Spectrum')
pl.xlabel('Wavelength (Angstroms)')
pl.ylabel('normalized flux / filter transmission')

#----------------------------------------------------------------------
# Plot filters in gray with several redshifted spectra
pl.figure()

redshifts = [0.0, 0.4, 0.8]
colors = 'bgr'

for z, c in zip(redshifts, colors):
    pl.plot((1. + z) * Xref[:, 0], Xref[:, 1], color=c)

pl.gca().add_patch(Arrow(4200, 0.47, 1300, 0, lw=0, width=0.05, color='r'))
pl.gca().add_patch(Arrow(5800, 0.47, 1250, 0, lw=0, width=0.05, color='r'))

pl.text(3800, 0.49, 'z = 0.0', fontsize=14, color=colors[0])
pl.text(5500, 0.49, 'z = 0.4', fontsize=14, color=colors[1])
pl.text(7300, 0.49, 'z = 0.8', fontsize=14, color=colors[2])

for f in 'ugriz':
    X = fetch_filter(f)
    pl.fill(X[:, 0], X[:, 1], ec='k', fc='k', alpha=0.2)

kwargs = dict(fontsize=20, color='gray', ha='center', va='center')
pl.text(3500, 0.02, 'u', **kwargs)
pl.text(4600, 0.02, 'g', **kwargs)
pl.text(6100, 0.02, 'r', **kwargs)
pl.text(7500, 0.02, 'i', **kwargs)
pl.text(8800, 0.02, 'z', **kwargs)

pl.xlim(3000, 11000)
pl.ylim(0, 0.55)

pl.title('Redshifting of a Spectrum')
pl.xlabel('Observed Wavelength (Angstroms)')
pl.ylabel('normalized flux / filter transmission')

pl.show()
