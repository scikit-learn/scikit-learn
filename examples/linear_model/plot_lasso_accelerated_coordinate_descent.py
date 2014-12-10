"""
=====================
Accelerated coordinate descent for the lasso
=====================

Lasso (L1  penalisation) implemented using an accelerated
coordinate descent, compared with coordinate descent.

"""
print(__doc__)

# Author: Olivier Fercoq <olivier.fercoq@telecom-paristech.fr>
# License: BSD 3 clause

import numpy as np
from scipy import linalg
from scipy import sparse
import matplotlib.pyplot as plt

import time

from sklearn.linear_model import lasso_path, enet_path

from sklearn.datasets.mldata import fetch_mldata
from sklearn.datasets.covtype import fetch_covtype
from sklearn.datasets.species_distributions import fetch_species_distributions

import os
from sklearn.externals import joblib, six
import zipfile

if six.PY3:
    from urllib.request import urlopen
else:
    from urllib2 import urlopen

from sklearn.utils.sparsefuncs import mean_variance_axis, inplace_column_scale

#data = fetch_mldata('housing')
#data = fetch_mldata('diabetes')
#data = fetch_mldata('regression-datasets-wisconsin')
#data = fetch_mldata('iris');
#data = fetch_covtype()

precompute = False

for experiment in range(4):
    if experiment == 0:
        data = fetch_mldata('leukemia')
        X = data.data
        y = data.target
        
        print("Comparing coordinate descent and accelerated"
              " coordinate descent on the leukemia dataset."
              " This is a dense dataset.")

    if experiment > 0:
        # Download BlogFeedback dataset
        URL = ("http://archive.ics.uci.edu/ml/machine-learning-databases"
               "/00304/BlogFeedback.zip")
        target_dir = "../../../../scikit_learn_data/BlogFeedback"
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
    
        archive_name = "BlogFeedback.zip"
        file_name = "blogData_train.csv"
        archive_path = os.path.join(target_dir, archive_name)
        file_path = os.path.join(target_dir, file_name)
        if not os.path.exists(file_path):
            print("Downloading dataset (2.5 MB)")
            opener = urlopen(URL)
            open(archive_path, 'wb').write(opener.read())

            zf = zipfile.ZipFile(archive_path)
            for i, name in enumerate(zf.namelist()):
                open(file_path, 'wb').write(zf.read(name))
                #we just keep the last one, which is the train dataset

        data = np.loadtxt(file_path, delimiter=',');
        X = data[:,:-1];
        y = data[:,X.shape[1]]
        
    if experiment == 1:
        print("Comparing coordinate descent and accelerated"
              " coordinate descent on the BlogFeedback dataset."
              " We precompute the Gram matrix.")
        precompute = True
    if experiment == 2:      
        print("Comparing coordinate descent and accelerated"
              " coordinate descent on the BlogFeedback dataset."
              " This is a dense dataset and we do not precompute"
              " the Gram matrix (long experiment).")
        precompute = False
    if experiment == 3:
        print("Comparing coordinate descent and accelerated"
              " coordinate descent on the BlogFeedback dataset."
              " We consider the data matrix as a sparse matrix (long experiment).")
        X = sparse.csc_matrix(X)
        
    X = X.astype(float)
    y = y.astype(float)

    # Standardize data (easier to set the l1_ratio parameter)
    y /= linalg.norm(y)
    if sparse.issparse(X) is not True:
        X_mean = X.mean(axis=0)
        X_std = X.std(axis=0)
        X = X - X_mean
        X /= np.maximum(1e-10, X_std) 
    else:
        X_mean, X_var = mean_variance_axis(X, axis=0)
        X_var *= X.shape[0]
        X_std = np.sqrt(X_var, X_var)
        X_std[X_std == 0] = 1
        inplace_column_scale(X, 1. / X_std)

    # Compute paths
    eps = 1e-3 # the smaller it is the longer is the path

    tols = range(4, 9, 4)
    accels = {True, False}
    times = np.empty((len(accels), len(tols)))


    plt.close('all')

    for itol, tol in enumerate(tols):
        plt.figure()
        for iaccel, accel in enumerate(accels):
            begin = time.time()
            print("Computing the regularization path using the lasso with tol = %g and acceleration = %d" %(10 ** (-tol), accel) )
            alphas_lasso, coefs_lasso, gap_lasso = lasso_path(X, y, eps, n_alphas=10,
                                 verbose=2, max_iter=20000/(1+accel), 
                                 accelerated=accel, tol=10 ** (-tol),
                                 selection='random', fit_intercept=False, normalize=False,
                                 X_mean=X_mean, X_std=X_std, precompute=precompute)
            # We allow accelerated coordinate descent half the number of iteration
            # because each iteration is twice as long

            duration = time.time() - begin
            print(duration)


            times[iaccel, itol] = duration

            # Did the algorithm converge ?
            gap_lasso = np.maximum(np.abs(gap_lasso), 1e-15)

            l1 = plt.plot(-np.log10(alphas_lasso), np.log10(gap_lasso.T), label=accels)
            plt.axhline(-tol, linestyle='--', color='k')

        plt.xlabel('-Log(alpha)')
        plt.ylabel('Log(gap)')
        plt.title('Did the algorithm converge on experiment %d?' %experiment)
        plt.axis('tight')
        plt.legend([ 'classical', 'goal','accelerated'])

    # Which algorithm is faster ?
    plt.figure()
    plt.title('Which algorithm is faster on experiment %d?' %experiment)
    plt.plot(tols, times.T)
    plt.legend(['classical', 'accelerated'])
    plt.xlabel("Accuracy")
    plt.ylabel("Time (s)")
    plt.show()

