"""
=============================
Ledoit-Wolf vs OAS estimation
=============================

"""
print __doc__

import numpy as np
import pylab as pl

from scikits.learn.covariance import LedoitWolf, OAS

###############################################################################
# Generate sample data
n_features = 100

n_samples_range = np.arange(6, 30, 1)
repeat = 50
lw_mse = np.zeros((n_samples_range.size, repeat))
oa_mse = np.zeros((n_samples_range.size, repeat))
lw_shrinkage = np.zeros((n_samples_range.size, repeat))
oa_shrinkage = np.zeros((n_samples_range.size, repeat))
for i, n_samples in enumerate(n_samples_range):
    for j in range(repeat):
        coloring_matrix = np.eye(n_features)
        X_train = np.dot(
            np.random.normal(size=(n_samples, n_features)), coloring_matrix)
        X_test = np.dot(
            np.random.normal(size=(n_samples, n_features)), coloring_matrix)
        
        lw = LedoitWolf(store_precision=False)
        lw.fit(X_train)
        lw_mse[i,j] = lw.mse(np.eye(n_features))
        lw_shrinkage[i,j] = lw.shrinkage_

        oa = OAS(store_precision=False)
        oa.fit(X_train)
        oa_mse[i,j] = oa.mse(np.eye(n_features))
        oa_shrinkage[i,j] = oa.shrinkage_
        
pl.subplot(2,1,1)
pl.errorbar(n_samples_range, lw_mse.mean(1), yerr=lw_mse.std(1),
            label='Ledoit-Wolf', color='g')
pl.errorbar(n_samples_range, oa_mse.mean(1), yerr=oa_mse.std(1),
            label='OAS', color='r')
pl.xlabel("n_samples")
pl.ylabel("MSE")
pl.legend(loc="upper right")
pl.title("Comparison of covariance estimators")

pl.subplot(2,1,2)
pl.errorbar(n_samples_range, lw_shrinkage.mean(1), yerr=lw_shrinkage.std(1),
            label='Ledoit-Wolf', color='g')
pl.errorbar(n_samples_range, oa_shrinkage.mean(1), yerr=oa_shrinkage.std(1),
            label='OAS', color='r')
pl.xlabel("n_samples")
pl.ylabel("Shrinkage")
pl.legend(loc="lower right")
pl.ylim(pl.ylim()[0], 1. + (pl.ylim()[1] - pl.ylim()[0])/10.)

pl.show()
