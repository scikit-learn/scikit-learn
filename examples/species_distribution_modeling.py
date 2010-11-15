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
density estimation problem and use the OneClassSVM provided
by the package `scikits.learn.svm` as our modeling tool.
The data that we use in this example comes from Phillips et. al. (2006).

The two species are:

 - `Bradypus variegatus <http://en.wikipedia.org/wiki/Bradypus_variegatus>`_ ,
   the Brown-throated Sloth.

 - `Microryzomys minutus <http://en.wikipedia.org/wiki/Microryzomys_minutus>`_ ,
   also known as the Forest Small Rice Rat, a rodent that lives in Peru,
   Colombia, Ecuador, Peru, and Venezuela.

Some of the environmental variables (=coverages) are 
listed below. Note that we do not use any spacial variables.

Environmental variables
-----------------------
Annual means from 1961 to 1990.
dtr6190_ann Diurnal temperature range
frs6190_ann Ground frost frequency
pre6190_ann Precipitation
tmn6190_ann Minimum temperature
tmp6190_ann Mean temperature
tmx6190_ann Maximum temperature
vap6190_ann Vapor pressure
frs6190_ann Ground frost frequency
pre6190_ann Precipitation

References
----------

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


def load_dir(directory):
    data = []
    for fpath in glob("%s/*.asc" % normpath(directory)):
        fname = split(fpath)[-1]
        fname = fname[:fname.index(".")]
        X = np.loadtxt(fpath, skiprows=6)
        data.append((fname, X))
    return dict(data)


def get_coverages(points, coverages, xx, yy):
    X = np.empty((points.shape[0], len(coverages)))
    for n in range(points.shape[0]):
        i = np.searchsorted(xx, points[n, 0])
        j = np.searchsorted(yy, points[n, 1])
        for d, fname in enumerate(coverages):
            X[n, d] = coverages[fname][-j, i]
    return X

################################################################################
# Download the data, if not already on disk
samples_url = "http://www.cs.princeton.edu/~schapire/maxent/datasets/samples.zip"
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
ncols = 1212
nrows = 1592
xllcorner = -94.8
yllcorner = -56.05
grd_sz = 0.05

# x,y coordinates for each cell
xmin = xllcorner + grd_sz
xmax = xmin + (ncols * grd_sz)
ymin = yllcorner + grd_sz
ymax = ymin + (nrows * grd_sz)
xx = np.arange(xmin, xmax, grd_sz)
yy = np.arange(ymin, ymax, grd_sz)

print "Data grid"
print "---------"
print "xmin, xmax:", xmin, xmax
print "ymin, ymax:", ymin, ymax
print "grid size:", grd_sz

# Per species data
train_bv = train[train[:,0] == 0][:, 1:]
test_bv = test[test[:,0] == 0][:, 1:]
train_mm = train[train[:,0] == 1][:, 1:]
test_mm = test[test[:,0] == 1][:, 1:]

# Load env variable grids
coverage = load_dir("coverages")

# Get features (=coverages)
print "compute features (=coverages)...",
train_bv_cover = get_coverages(train_bv, coverage, xx, yy)
test_bv_cover = get_coverages(test_bv, coverage, xx, yy)
train_mm_cover = get_coverages(train_mm, coverage, xx, yy)
test_mm_cover = get_coverages(test_mm, coverage, xx, yy)
print "done."

# Compute coverage for whole grid.
X, Y = np.meshgrid(xx, yy)
points = np.c_[X.ravel(), Y.ravel()]
points_data = get_coverages(points, coverage, xx, yy)
#basemap = False
for i, species, cover, observations in zip(count(), species,
                                           [train_bv_cover, train_mm_cover],
                                           [train_bv, train_mm]):
    print "_" * 80
    print "Modeling distribution of species '%s'" % species
    print
    # Standardize features
    mean = cover.mean(axis=0)
    std = cover.std(axis=0)
    cover = (cover - mean) / std

    # Fit OneClassSVM
    print "fit OneCLassSVM ... ", 
    clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.5)
    clf.fit(cover)
    print "done. "

    # Plot map of South America
    pl.subplot(1,2,i + 1)
    if basemap:
        print "plot coastlines using basemap"
        m = Basemap(projection='cyl', llcrnrlat=ymin,
                urcrnrlat=ymax, llcrnrlon=xmin,
                urcrnrlon=xmax, resolution='c')
        m.drawcoastlines()
        m.drawcountries()
        m.drawlsmask(land_color='0.8', ocean_color='w')
    else:
        print "plot coastlines from cover"
        tmp = coverage['tmp6190_ann']
        XX, YY = np.meshgrid(xx, yy[::-1])
        Z = np.empty(XX.shape)
        for (row, col), val in np.ndenumerate(XX):
            Z[row, col] = tmp[row, col]
        CS = pl.contour(XX, YY, Z, levels=[-9999])

    points_data_std = (points_data - mean) / std
    print "predicting points..."
    Z = clf.predict_margin(points_data_std)
    levels = np.linspace(Z.min(), Z.max(), 25)
    # Blank out grid cells on the sea.
    Z[points_data[:,2] == -9999] = -9999
    Z = Z.reshape(X.shape)
    CS = pl.contourf(X, Y, Z, levels=levels, cmap=pl.cm.Reds)
    pl.colorbar()
    pl.scatter(observations[:, 0], observations[:, 1], s=8, c='black',
               marker='^', label='observations')
    pl.legend()
    pl.title(species)

pl.show()
