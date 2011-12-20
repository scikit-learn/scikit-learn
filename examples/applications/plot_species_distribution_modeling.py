"""
=============================
Species distribution modeling
=============================

Modeling species' geographic distributions is an important
problem in conservation biology. In this example we
model the geographic distribution of two south american
mammals given past observations and 14 environmental
variables. Since we have only positive examples (there are
no unsuccessful observations), we cast this problem as a
density estimation problem and use the `OneClassSVM` provided
by the package `sklearn.svm` as our modeling tool.
The dataset is provided by Phillips et. al. (2006).
If available, the example uses
`basemap <http://matplotlib.sourceforge.net/basemap/doc/html/>`_
to plot the coast lines and national boundaries of South America.

The two species are:

 - `"Bradypus variegatus"
   <http://www.iucnredlist.org/apps/redlist/details/3038/0>`_ ,
   the Brown-throated Sloth.

 - `"Microryzomys minutus"
   <http://www.iucnredlist.org/apps/redlist/details/13408/0>`_ ,
   also known as the Forest Small Rice Rat, a rodent that lives in Peru,
   Colombia, Ecuador, Peru, and Venezuela.

References:

 * `"Maximum entropy modeling of species geographic distributions"
   <http://www.cs.princeton.edu/~schapire/papers/ecolmod.pdf>`_
   S. J. Phillips, R. P. Anderson, R. E. Schapire - Ecological Modelling,
   190:231-259, 2006.

"""

from cStringIO import StringIO
import os
from time import time

import numpy as np
import pylab as pl
from scipy.sparse import csr_matrix

from sklearn.externals import joblib
from sklearn.datasets.base import Bunch
from sklearn.datasets.species_distributions import fetch_species_distributions
from sklearn import svm
from sklearn.metrics import roc_curve, auc

# if basemap is available, we'll use it.
# otherwise, we'll improvise later...
try:
    from mpl_toolkits.basemap import Basemap
    basemap = True
except ImportError:
    basemap = False


def create_species_bunch(species_name,
                         train, test,
                         coverages, xgrid, ygrid):
    """
    create a bunch with information about a particular organism

    This will use the test/train record arrays to extract the
    data specific to the given species name.
    """
    bunch = Bunch(name=' '.join(species_name.split("_")[:2]))

    points = dict(test=test, train=train)

    for label, pts in points.iteritems():
        # choose points associated with the desired species
        pts = pts[pts['species'] == species_name]
        bunch['pts_%s' % label] = pts

        # determine coverage values for each of the training & testing points
        ix = np.searchsorted(xgrid, pts['dd long'])
        iy = np.searchsorted(ygrid, pts['dd lat'])
        bunch['cov_%s' % label] = coverages[:, -iy, ix].T

    return bunch
            

def plot_species_distribution(species = ["bradypus_variegatus_0",
                                         "microryzomys_minutus_0"]):
    """
    Plot the species distribution.
    """
    t0 = time()

    # Load the compressed data
    B = fetch_species_distributions()

    # Set up the data
    species_map = dict([(s, i) for i, s in enumerate(species)])

    # x,y coordinates for corner cells
    xmin = B.x_left_lower_corner + B.grid_size
    xmax = xmin + (B.Nx * B.grid_size)
    ymin = B.y_left_lower_corner + B.grid_size
    ymax = ymin + (B.Ny * B.grid_size)

    # x coordinates of the grid cells
    xgrid = np.arange(xmin, xmax, B.grid_size)
    # y coordinates of the grid cells
    ygrid = np.arange(ymin, ymax, B.grid_size)

    # The grid in x,y coordinates
    X, Y = np.meshgrid(xgrid, ygrid[::-1])
    
    # convert coverages to dense array
    #coverages = np.asarray([mat.toarray() for mat in B.coverages],
    #                       dtype=np.float32)
    coverages = B.coverages
    
    # create a bunch for each species
    BV = create_species_bunch(species[0],
                              B.train, B.test,
                              coverages, xgrid, ygrid)
    MM = create_species_bunch(species[1],
                              B.train, B.test,
                              coverages, xgrid, ygrid)

    # background points (grid coordinates) for evaluation
    np.random.seed(13)
    background_points = np.c_[np.random.randint(low=0, high=B.Ny,
                                                size=10000),
                              np.random.randint(low=0, high=B.Nx,
                                                size=10000)].T

    # Fit, predict, and plot for each species.
    for i, species in enumerate([BV, MM]):
        print "_" * 80
        print "Modeling distribution of species '%s'" % species.name
        print
        # Standardize features
        mean = species.cov_train.mean(axis=0)
        std = species.cov_train.std(axis=0)
        train_cover_std = (species.cov_train - mean) / std

        # Fit OneClassSVM
        print "fit OneClassSVM ... ",
        clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.5)
        clf.fit(train_cover_std)
        print "done. "

        # Plot map of South America
        pl.subplot(121 + i)
        if basemap:
            print "plot coastlines using basemap"
            m = Basemap(projection='cyl', llcrnrlat=ymin,
                        urcrnrlat=ymax, llcrnrlon=xmin,
                        urcrnrlon=xmax, resolution='c')
            m.drawcoastlines()
            m.drawcountries()
            m.drawrivers()
        else:
            print "plot coastlines from coverage"
            CS = pl.contour(X, Y, coverages[2], levels=[-9999], colors="k",
                            linestyles="solid")
            pl.xticks([])
            pl.yticks([])
        
        print "predict species distribution"
        
        # Predict species distribution using the training data
        Z = np.ones((B.Ny, B.Nx), dtype=np.float64)
        
        # find the land points
        idx = np.where(coverages[2] > -9999)

        coverages_land = coverages[:, idx[0], idx[1]].T

        pred = clf.decision_function((coverages_land - mean) / std)[:, 0]
        Z *= pred.min()
        Z[idx[0], idx[1]] = pred

        levels = np.linspace(Z.min(), Z.max(), 25)
        Z[coverages[2] == -9999] = -9999

        # plot contours of the prediction
        CS = pl.contourf(X, Y, Z, levels=levels, cmap=pl.cm.Reds)
        pl.colorbar(format='%.2f')

        # scatter training/testing points
        pl.scatter(species.pts_train['dd long'], species.pts_train['dd lat'],
                   s=2 ** 2, c='black',
                   marker='^', label='train')
        pl.scatter(species.pts_test['dd long'], species.pts_test['dd lat'],
                   s=2 ** 2, c='black',
                   marker='x', label='test')
        pl.legend()
        pl.title(species.name)
        pl.axis('equal')

        # Compute AUC w.r.t. background points
        pred_background = Z[background_points[0], background_points[1]]
        pred_test = clf.decision_function((species.cov_test - mean) / std)[:, 0]
        scores = np.r_[pred_test, pred_background]
        y = np.r_[np.ones(pred_test.shape), np.zeros(pred_background.shape)]
        fpr, tpr, thresholds = roc_curve(y, scores)
        roc_auc = auc(fpr, tpr)
        pl.text(-35, -70, "AUC: %.3f" % roc_auc, ha="right")
        print "Area under the ROC curve : %f" % roc_auc

    print "time elapsed: %.3fs" % (time() - t0)

    pl.show()

if __name__ == '__main__':
    plot_species_distribution()
    pl.show()
