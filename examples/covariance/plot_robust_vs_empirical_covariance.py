"""
=======================================
Robust vs Empirical covariance estimate
=======================================

The Minimum Covariance Determinant estimator (MCD) is a robust
estimator of a data set's covariance introduced by P.J.Rousseuw in
[1]. The idea is to find a given proportion of "good" observations
which are not outliers and compute their empirical covariance matrix.
From these observations, it is necessary to compute a robust location
estimate before computing the Minimal Covariance Determinant estimator.

In this example, we compare the estimation errors that are made when
using four types of location and covariance estimates on contaminated
gaussian distributed data sets:
 - The mean and the empirical covariance of the full dataset, which break
   down as soon as there are outliers in the data set
 - The robust MCD, that has a low error provided n_samples > 5 * n_features
 - The robust MCD reweighted according to Rousseeuw's recommandations:
   observations are given a 0 weight if they are found to be outlying
   according to their MCD-based Mahalanobis distance. Doing so improve the
   efficiency of the estimators at gaussian models.
 - The mean and the empirical covariance of the observations that are known
   to be good ones. This can be considered as a "perfect" MCD estimation.

[1] P. J. Rousseeuw. Least median of squares regression. J. Am
    Stat Ass, 79:871, 1984.

"""
print __doc__

import numpy as np
import pylab as pl
import matplotlib.font_manager

from scikits.learn.covariance import EmpiricalCovariance, MCD

n_samples = 100
n_features = 5
repeat = 10

range_n_outliers = np.arange(0, n_samples / 2, 5)
error_location_mcd = np.zeros((range_n_outliers.size, repeat))
error_covariance_mcd = np.zeros((range_n_outliers.size, repeat))
error_location_mcdr = np.zeros((range_n_outliers.size, repeat))
error_covariance_mcdr = np.zeros((range_n_outliers.size, repeat))
error_location_emp_full = np.zeros((range_n_outliers.size, repeat))
error_covariance_emp_full = np.zeros((range_n_outliers.size, repeat))
error_location_emp_pure = np.zeros((range_n_outliers.size, repeat))
error_covariance_emp_pure = np.zeros((range_n_outliers.size, repeat))
for i, n_outliers in enumerate(range_n_outliers):
    for j in range(repeat):
        # generate data
        data = np.random.randn(n_samples, n_features)
        # add some outliers
        outliers_index = np.random.permutation(n_samples)[:n_outliers]
        outliers_offset = 10. * \
            (np.random.randint(2, size=(n_outliers,n_features)) - 0.5)
        data[outliers_index] += outliers_offset
        inliers_mask = np.ones(n_samples).astype(bool)
        inliers_mask[outliers_index] = False
        
        # fit a Minimum Covariance Determinant (MCD) robust estimator to data
        S = MCD().fit(data, reweight=None)
        # compare robust estimates with the true location and covariance
        error_location_mcd[i,j] = np.sum(S.location_ ** 2)
        error_covariance_mcd[i,j] = S.error(
            np.eye(n_features), error_type='rmse')
        # fit a reweighted MCD robust estimator to data
        S = MCD().fit(data)
        # compare robust estimates with the true location and covariance
        error_location_mcdr[i,j] = np.sum(S.location_ ** 2)
        error_covariance_mcdr[i,j] = S.error(
            np.eye(n_features), error_type='rmse')
        # compare estimators learnt from the full data set with true parameters
        error_location_emp_full[i,j] = np.sum(data.mean(0) ** 2)
        error_covariance_emp_full[i,j] = EmpiricalCovariance().fit(data).error(
            np.eye(n_features), error_type='rmse')
        # compare with an empirical covariance learnt from a pure data set
        # (i.e. "perfect" MCD)
        pure_data = data[inliers_mask]
        pure_location = pure_data.mean(0)
        pure_emp_cov = EmpiricalCovariance().fit(pure_data)
        error_location_emp_pure[i,j] = np.sum(pure_location ** 2)
        error_covariance_emp_pure[i,j] = pure_emp_cov.error(
            np.eye(n_features), error_type='rmse')    

# Display results
font_prop = matplotlib.font_manager.FontProperties(size=12)
pl.subplot(2,1,1)
pl.errorbar(range_n_outliers, error_location_mcd.mean(1),
            yerr=error_location_mcd.std(1)/np.sqrt(repeat),
            label="Robust location", color='cyan')
pl.errorbar(range_n_outliers, error_location_mcdr.mean(1),
            yerr=error_location_mcdr.std(1)/np.sqrt(repeat),
            label="Robust location (reweighted)", color='blue')
pl.errorbar(range_n_outliers, error_location_emp_full.mean(1),
            yerr=error_location_emp_full.std(1)/np.sqrt(repeat),
            label="Full data set mean", color='green')
pl.errorbar(range_n_outliers, error_location_emp_pure.mean(1),
            yerr=error_location_emp_pure.std(1)/np.sqrt(repeat),
            label="Pure data set mean", color='black')
pl.title("Influence of outliers on the location estimation")
pl.ylabel(r"Error ($||\mu - \hat{\mu}||_2^2$)")
pl.legend(loc="upper left", prop=font_prop)

pl.subplot(2,1,2)
x_size = range_n_outliers.size
pl.errorbar(range_n_outliers, error_covariance_mcd.mean(1),
            yerr=error_covariance_mcd.std(1),
            label="Robust covariance (MCD)", color='cyan')
pl.errorbar(range_n_outliers, error_covariance_mcdr.mean(1),
            yerr=error_covariance_mcdr.std(1),
            label="Robust covariance (reweighted MCD)", color='blue')
pl.errorbar(range_n_outliers[:(x_size/4 + 1)],
            error_covariance_emp_full.mean(1)[:(x_size/4 + 1)],
            yerr=error_covariance_emp_full.std(1)[:(x_size/4 + 1)],
            label="Full data set empirical covariance", color='green')
pl.plot(range_n_outliers[(x_size/4):(x_size/2 - 1)],
         error_covariance_emp_full.mean(1)[(x_size/4):(x_size/2 - 1)],
         color='green', ls='--')
pl.errorbar(range_n_outliers, error_covariance_emp_pure.mean(1),
            yerr=error_covariance_emp_pure.std(1),
            label="Pure data set empirical covariance", color='black')
pl.title("Influence of outliers on the covariance estimation")
pl.xlabel("Amount of contamination (%)")
pl.ylabel("RMSE") 
pl.legend(loc="upper right", prop=font_prop)

pl.show()
