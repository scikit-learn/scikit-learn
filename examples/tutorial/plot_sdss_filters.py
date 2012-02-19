"""
Plot the SDSS filter curves over a reference spectrum
"""
import os
import urllib2

import numpy as np
import pylab as pl

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
pl.plot(Xref[:, 0], Xref[:, 1], '-k', label='Vega spectrum', lw=2)

for f,c in zip('ugriz', 'bgrmk'):
    X = fetch_filter(f)
    pl.fill(X[:, 0], X[:, 1], ec=c, fc=c, alpha=0.4, label='%s-band' % f)

pl.xlim(3000, 11000)
pl.legend()

pl.title('SDSS Filters and Reference Spectrum')
pl.xlabel('Wavelength (Angstroms)')
pl.ylabel('flux (normalized)')
pl.show()
