import os
import urllib
import numpy as np
import pylab as pl

from sklearn.datasets import get_data_home
from sklearn.neighbors import KNeighborsRegressor

URL = 'http://cas.sdss.org/public/en/tools/search/x_sql.asp'
ARCHIVE_FILE = 'sdss_galaxy_colors.npy'

N_objects = 50000

DTYPE = [('u', float),
         ('g', float),
         ('r', float),
         ('i', float),
         ('z', float),
         ('specClass', int),
         ('redshift', float),
         ('redshift_err', float)]

QUERY_TEXT = ('\n'.join(
        ("SELECT TOP %i" % N_objects,
         "   p.u, p.g, p.r, p.i, p.z, s.specClass, s.z, s.zerr",
         "FROM PhotoObj AS p",
         "   JOIN SpecObj AS s ON s.bestobjid = p.objid",
         "WHERE ",
         "   p.u BETWEEN 0 AND 19.6",
         "   AND p.g BETWEEN 0 AND 20",
         "   AND s.specClass > 1 -- not UNKNOWN or STAR",
         "   AND s.specClass <> 5 -- not SKY",
         "   AND s.specClass <> 6 -- not STAR_LATE")))

def sql_query(sql_str, url=URL, format='csv'):
    """Execute SQL query"""
    # remove comments from string
    sql_str = ' \n'.join(map(lambda x: x.split('--')[0], sql_str.split('\n')))
    params = urllib.urlencode(dict(cmd=sql_str, format=format))
    return urllib.urlopen(url + '?%s' % params)


def fetch_photoz_data(data_home=None):
    data_home = get_data_home(data_home)

    archive_file = os.path.join(data_home, ARCHIVE_FILE)

    if not os.path.exists(archive_file):
        print "querying for %i objects" % N_objects
        print QUERY_TEXT
        output = sql_query(QUERY_TEXT)
        print "finished.  Processing & saving data"

        data = np.loadtxt(output, delimiter=',', skiprows=1, dtype=DTYPE)
        np.save(archive_file, data)
    else:
        data = np.load(archive_file)

    return data


data = fetch_photoz_data()

N = len(data)

# shuffle data
np.random.seed(0)
np.random.shuffle(data)

# put colors in a matrix
X = np.zeros((N, 4))
X[:, 0] = data['u'] - data['g']
X[:, 1] = data['g'] - data['r']
X[:, 2] = data['r'] - data['i']
X[:, 3] = data['i'] - data['z']
z = data['redshift']

# divide into training and testing data
Ntrain = N/2
Xtrain = X[:Ntrain]
ztrain = z[:Ntrain]

Xtest = X[Ntrain:]
ztest = z[Ntrain:]

knn = KNeighborsRegressor(1, weights='uniform')
zpred = knn.fit(Xtrain, ztrain).predict(Xtest)

axis_lim = np.array([-0.1, 2.5])

rms = np.sqrt(np.mean((ztest - zpred) ** 2))
print rms
print len(ztest)
print np.sum(abs(ztest - zpred) > 1)

ax = pl.axes()
pl.scatter(ztest, zpred, c='k', lw=0, s=4)
pl.plot(axis_lim, axis_lim, '--k')
pl.plot(axis_lim, axis_lim + rms, ':k')
pl.plot(axis_lim, axis_lim - rms, ':k')
pl.xlim(axis_lim)
pl.ylim(axis_lim)

pl.text(0.99, 0.02, "RMS error = %.2g" % rms,
        ha='right', va='bottom', transform=ax.transAxes,
        bbox=dict(ec='w', fc='w'), fontsize=16)

pl.title('Photo-z: Nearest Neigbor Regression')
pl.xlabel(r'$\mathrm{z_{true}}$', fontsize=14)
pl.ylabel(r'$\mathrm{z_{phot}}$', fontsize=14)
pl.show()



