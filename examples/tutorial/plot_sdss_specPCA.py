"""
SDSS Spectra Plots
------------------

This plots some of the SDSS spectra examples for the astronomy tutorial
"""
import os
import urllib2

import numpy as np
import pylab as pl

from sklearn.datasets import get_data_home
from sklearn import preprocessing
from sklearn.decomposition import RandomizedPCA

DATA_URL = ('http://www.astro.washington.edu/users/'
            'vanderplas/pydata/spec4000_corrected.npz')

def fetch_sdss_spec_data(data_home=None):
    data_home = get_data_home(data_home)

    local_file = os.path.join(data_home, os.path.basename(DATA_URL))

    # data directory is password protected so the public can't access it    
    password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
    password_mgr.add_password(None, DATA_URL, 'pydata', 'astroML')
    handler = urllib2.HTTPBasicAuthHandler(password_mgr)
    opener = urllib2.build_opener(handler)

    # download training data
    if not os.path.exists(local_file):
        fhandle = opener.open(DATA_URL)
        open(local_file, 'w').write(fhandle.read())

    return np.load(local_file)

#----------------------------------------------------------------------
#
#  Load the data
data = fetch_sdss_spec_data()

wavelengths = data['wavelengths']
X = data['X']
y = data['y']
labels = data['labels']

from matplotlib.ticker import FuncFormatter
format = FuncFormatter(lambda i, *args: labels[i].replace(' ', '\n'))

#----------------------------------------------------------------------
#
#  Plot the first few spectra, offset so they don't overlap
#
pl.figure()

for i_class in (2, 3, 4, 5, 6):
    i = np.where(y == i_class)[0][0]
    l = pl.plot(wavelengths, X[i] + 20 * i_class)
    c = l[0].get_color()
    pl.text(6800, 2 + 20 * i_class, labels[i_class], color=c)
            
pl.subplots_adjust(hspace=0)
pl.xlabel('wavelength (Angstroms)')
pl.ylabel('flux + offset')
pl.title('Sample of Spectra')

#----------------------------------------------------------------------
#
#  Plot the mean spectrum
#
X = preprocessing.normalize(X, 'l2')

pl.figure()

mu = X.mean(0)
std = X.std(0)

pl.plot(wavelengths, mu, color='black')
pl.fill_between(wavelengths, mu - std, mu + std, color='#CCCCCC')
pl.xlim(wavelengths[0], wavelengths[-1])
pl.ylim(0, 0.06)
pl.xlabel('wavelength (Angstroms)')
pl.ylabel('scaled flux')
pl.title('Mean Spectrum')

#----------------------------------------------------------------------
#
#  Plot a random pair of digits
#
pl.figure()
np.random.seed(25255)
i1, i2 = np.random.randint(1000, size=2)

pl.scatter(X[:, i1], X[:, i2], c=y, s=4, lw=0,
           vmin=2, vmax=6, cmap=pl.cm.jet)
pl.colorbar(ticks = range(2, 7), format=format)
pl.xlabel('wavelength = %.1f' % wavelengths[i1])
pl.ylabel('wavelength = %.1f' % wavelengths[i2])
pl.title('Random Pair of Spectra Bins')

#----------------------------------------------------------------------
#
#  Perform PCA
#

rpca = RandomizedPCA(n_components=4, random_state=0)
X_proj = rpca.fit_transform(X)

#----------------------------------------------------------------------
#
#  Plot PCA components
#

pl.figure()
pl.scatter(X_proj[:, 0], X_proj[:, 1], c=y, s=4, lw=0,
           vmin=2, vmax=6, cmap=pl.cm.jet)
pl.colorbar(ticks = range(2, 7), format=format)
pl.xlabel('coefficient 1')
pl.ylabel('coefficient 2')
pl.title('PCA projection of Spectra')

#----------------------------------------------------------------------
#
#  Plot PCA eigenspectra
#

pl.figure()

l = pl.plot(wavelengths, rpca.mean_ - 0.15)
c = l[0].get_color()
pl.text(7000, -0.16, "mean" % i, color=c)

for i in range(4):
    l = pl.plot(wavelengths, rpca.components_[i] + 0.15 * i)
    c = l[0].get_color()
    pl.text(7000, -0.01 + 0.15 * i, "component %i" % (i + 1), color=c)
pl.ylim(-0.2, 0.6)
pl.xlabel('wavelength (Angstroms)')
pl.ylabel('scaled flux + offset')
pl.title('Mean Spectrum and Eigen-spectra')

pl.show()
