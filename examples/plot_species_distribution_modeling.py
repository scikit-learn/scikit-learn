"""
=============================
Species distribution modeling
=============================

Modeling species geographic distributions is an important
problem in conservation biology. In this example we
model the geographic distribution of two south american
mammals given past observations and 14 environmental
variables. Since we have only positive examples (there are
no unsuccessful observations), we cast this problem as a
density estimation problem and use the `OneClassSVM` provided
by the package `scikits.learn.svm` as our modeling tool.
The data that we use in this example comes from Phillips et. al. (2006).
If available, the example uses `basemap <http://matplotlib.sourceforge.net/basemap/doc/html/>`_
to plot the coast lines and national boundaries of south america.

The two species are:

 - `Bradypus variegatus <http://www.iucnredlist.org/apps/redlist/details/3038/0>`_ ,
   the Brown-throated Sloth.

 - `Microryzomys minutus <http://www.iucnredlist.org/apps/redlist/details/13408/0>`_ ,
   also known as the Forest Small Rice Rat, a rodent that lives in Peru,
   Colombia, Ecuador, Peru, and Venezuela.
   

Some of the environmental variables (=coverages) are
listed below. Note that we do not use any spacial variables.

Environmental variables:

Annual means from 1961 to 1990.
 - dtr6190_ann Diurnal temperature range
 - frs6190_ann Ground frost frequency
 - pre6190_ann Precipitation
 - tmn6190_ann Minimum temperature
 - tmp6190_ann Mean temperature
 - tmx6190_ann Maximum temperature
 - vap6190_ann Vapor pressure
 - frs6190_ann Ground frost frequency
 - pre6190_ann Precipitation

References:

 * `"Maximum entropy modeling of species geographic distributions"
   <http://www.cs.princeton.edu/~schapire/papers/ecolmod.pdf>`_
   S. J. Phillips, R. P. Anderson, R. E. Schapire - Ecological Modelling,
   190:231-259, 2006.

"""
from __future__ import division

print __doc__

import pylab as pl
import numpy as np

try:
    from mpl_toolkits.basemap import Basemap
    basemap = True
except ImportError:
    basemap = False

from os.path import normpath, split, exists
from glob import glob
from itertools import count
from scikits.learn import svm
from time import time


def load_dir(directory):
    """Loads each of the grids and returns a
    tensor of shape [14, n_rows, n_cols].
    """
    data = []
    for fpath in glob("%s/*.asc" % normpath(directory)):
        fname = split(fpath)[-1]
        fname = fname[:fname.index(".")]
        X = np.loadtxt(fpath, skiprows=6)
        data.append(X)
    return np.array(data)


def get_coverages(points, coverages, xx, yy):
    """
    Returns
    -------
    array : shape = [n_points, 14]
    """
    rows = []
    cols = []
    for n in range(points.shape[0]):
        i = np.searchsorted(xx, points[n, 0])
        j = np.searchsorted(yy, points[n, 1])
        rows.append(-j)
        cols.append(i)
    return coverages[:, rows, cols].T

################################################################################
# Download the data, if not already on disk
samples_url = "http://www.cs.princeton.edu/~schapire/maxent/datasets/" \
              "samples.zip"
coverage_url = "http://www.cs.princeton.edu/~schapire/maxent/datasets/" \
               "coverages.zip"
samples_archive_name = "samples.zip"
coverage_archive_name = "coverages.zip"


def download(url, archive_name):
    if not exists(archive_name[:-4]):
        if not exists(archive_name):
            import urllib
            print "Downloading data, please wait ..."
            print url
            opener = urllib.urlopen(url)
            open(archive_name, 'wb').write(opener.read())
            print

        import zipfile
        print "Decompressiong the archive: " + archive_name
        zipfile.ZipFile(archive_name).extractall()
        print

t0 = time()

download(samples_url, samples_archive_name)
download(coverage_url, coverage_archive_name)

species = ["bradypus variegatus", "microryzomys minutus"]
species_map = dict([(s, i) for i, s in enumerate(species)])
species2id = lambda s: species_map[" ".join(s.split("_")[:2])]
train = np.loadtxt('samples/alltrain.csv', converters={0: species2id},
                   skiprows=1, delimiter=",")
test = np.loadtxt('samples/alltest.csv', converters={0: species2id},
                  skiprows=1, delimiter=",")

# Data resolution
n_cols = 1212
n_rows = 1592
x_left_lower_corner = -94.8
y_left_lower_corner = -56.05
grid_size = 0.05

# x,y coordinates for each cell
xmin = x_left_lower_corner + grid_size
xmax = xmin + (n_cols * grid_size)
ymin = y_left_lower_corner + grid_size
ymax = ymin + (n_rows * grid_size)
xx = np.arange(xmin, xmax, grid_size)
yy = np.arange(ymin, ymax, grid_size)

print "Data grid"
print "---------"
print "xmin, xmax:", xmin, xmax
print "ymin, ymax:", ymin, ymax
print "grid size:", grid_size

# Per species data
train_bv = train[train[:,0] == 0, 1:]
test_bv = test[test[:,0] == 0, 1:]
train_mm = train[train[:,0] == 1, 1:]
test_mm = test[test[:,0] == 1, 1:]

# Load env variable grids
coverage = load_dir("coverages")

# Get features (=coverages)
print "compute features (=coverages)...",
train_bv_cover = get_coverages(train_bv, coverage, xx, yy)
test_bv_cover = get_coverages(test_bv, coverage, xx, yy)
train_mm_cover = get_coverages(train_mm, coverage, xx, yy)
test_mm_cover = get_coverages(test_mm, coverage, xx, yy)
print "done."


def predict(clf, mean, std):
    """Predict the density of the land grid cells
    under the model `clf`.

    Returns
    -------
    array : shape [n_rows, n_cols]
    """
    Z = np.zeros((n_rows, n_cols))
    for row, col in np.ndindex(n_rows, n_cols):
        if coverage[2, row, col] > -9999:
            Z[row, col]= clf.predict_margin((coverage[:,row,col]-mean)/std)
    return Z

X, Y = np.meshgrid(xx, yy[::-1])
basemap = False
for i, species, cover, observations in zip(count(), species,
                                           [train_bv_cover, train_mm_cover],
                                           [train_bv, train_mm]):
    print "_" * 80
    print "Modeling distribution of species '%s'" % species
    print
    # Standardize features
    mean = cover.mean(axis=0)
    std = cover.std(axis=0)
    cover_std = (cover - mean) / std

    # Fit OneClassSVM
    print "fit OneClassSVM ... ",
    clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.5)
    clf.fit(cover_std)
    print "done. "

    # Plot map of South America
    pl.subplot(1, 2, i + 1)
    if basemap:
        print "plot coastlines using basemap"
        m = Basemap(projection='cyl', llcrnrlat=ymin,
                urcrnrlat=ymax, llcrnrlon=xmin,
                urcrnrlon=xmax, resolution='c')
        m.drawcoastlines()
        m.drawcountries()
        #m.drawrivers()
    else:
        print "plot coastlines from cover"
        CS = pl.contour(X, Y, coverage[2,:,:], levels=[-9999], colors="k", linestyles="solid")
        pl.xticks([])
        pl.yticks([])

    Z = predict(clf, mean, std)
    levels = np.linspace(Z.min(), Z.max(), 25)
    Z[coverage[2,:,:] == -9999] = Z.min()
    CS = pl.contourf(X, Y, Z, levels=levels, cmap=pl.cm.Reds)
    pl.colorbar(format='%.2f')
    pl.scatter(observations[:, 0], observations[:, 1], s=8, c='black',
               marker='^', label='observations')
    pl.legend()
    pl.title(species)
    pl.axis('equal')

print "time elapsed: %.3fs" % (time() - t0)

pl.show()
