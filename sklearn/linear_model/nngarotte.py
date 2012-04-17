"""
Non-Negative Garotte implementation with the Scikit-learn
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Jaques Grobler <jaques.grobler@inria.fr>
#
# License: BSD Style.


#imports
import numpy as np


from sklearn.linear_model.base import LinearModel
from sklearn.linear_model import LinearRegression, Lasso, lasso_path

#global functions--------------------------------------------------------
def non_negative_garotte(X, y, alpha, tol=0.001):
    """
    TODO - non_negative_garotte docstring
    """
    # Obtain the ordinary least squares coefficients from our data
    coef_ols = LinearRegression(fit_intercept=False).fit(X, y).coef_

    X = X * coef_ols[np.newaxis, :]
    # find the shrinkage factor by minimising the sum of square residuals
    # under the restriction that it is positive (positive=True)
    shrink_coef = Lasso(alpha=alpha, fit_intercept=False,
                        positive=True, normalize=False,
                        tol=tol).fit(X, y).coef_

    # Shrunken betas
    coef = coef_ols * shrink_coef

    # Residual Sum of Squares
    rss = np.sum((y - np.dot(X, coef)) ** 2)
    return coef, shrink_coef, rss


def non_negative_garotte_path(X, y, alpha):
    """
    TODO - non_negative_garotte_path docstring
    Compute the Non-negative Garotte path

    """

    # Obtain the ordinary least squares coefficients from our data
    coef_ols = LinearRegression(fit_intercept=False).fit(X, y).coef_

    X = X * coef_ols[np.newaxis, :]
    # find the shrinkage factor by minimising the sum of square residuals
    # under the restriction that it is positive (positive=True)
    # lasso_path returns a list of models - below is a bit of a hack
    # to get the coefficients of a model (all are identical if you fix
    # alpha.. Is there a better way to do this?
    shrink_coef = lasso_path(X, y, positive=True, alpha=alpha)[0].coef_
    
    # Shrunken betas
    coef_path = coef_ols * shrink_coef

    # Residual Sum of Squares
    rss = np.sum((y - np.dot(X, coef)) ** 2)
    
    return coef_path, shrink_coef, rss

#classes---------------------------------------------------------------------

class NonNegativeGarrote(LinearModel):
    """NonNegativeGarrote - TODO description

    Ref:
    Breiman, L. (1995), "Better Subset Regression Using the Nonnegative
    Garrote," Technometrics, 37, 373-384. [349,351]

    Parameters
    ----------
    TODO

    Attributes
    ----------
    TODO

    Examples
    --------
    TODO

    See also
    --------
    TODO



    NOTES:
    good default value for alpha?
    """
    def __init__(self, alpha=0.35, fit_intercept=True, tol=1e-4, normalize=False,
                 copy_X=True):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.tol = tol
        self.normalize = normalize
        self.copy_X = copy_X

    def fit(self, X, y):
        """Fit the model using X, y as training data.

        parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            training data.

        y : array-like, shape = [n_samples]
            target values.

        returns
        -------
        self : object
            returns an instance of self.
        """
        X, y, X_mean, y_mean, X_std = LinearModel._center_data(X, y,
                self.fit_intercept, self.normalize, self.copy_X)

        self.coef_, self.shrink_coef_, self.rss_ = \
                                    non_negative_garotte(X, y, alpha)
        self._set_intercept(X_mean, y_mean, X_std)



# For now, this is just a script for comparing the paths of the nng and the lasso
# with a synthetic dataset. Will make a proper example of this later
# This obtains the same figures as `figure 1` in 
# http://www2.isye.gatech.edu/statistics/papers/05-25.pdf , page 11

if __name__ == '__main__':
    import pylab as pl
    from sklearn.utils import check_random_state
    from sklearn.linear_model import lars_path
    

    print 'running nngarrote.py - __main__'
    rng = check_random_state(None)

    ng_path_correct = 0
    lars_path_correct = 0
    max_samples = 501
    
    # true path
    coef = np.array([1, 1, 0])
    
    # lists for plotting the two techniques results
    hits_lars = []
    hits_ng = []
    
    # for 4 different values of alpha
    for alpha_val, fig_num in ((0.35, 1), (0.45, 2), (0.55, 3), (0.65, 4)):
        print 'for alpha = ', alpha_val
        # set up plots
        pl.figure(fig_num, figsize=(5, 5))
        pl.clf
        pl.axis('tight')
        pl.title('alpha = %.2f' % alpha_val )
        pl.xlabel('Sample Size')
        pl.ylabel('Frequency of Selecting Correct Models')
        
        # vary the sample size from 25 up until 500 
        for sample_size in xrange(25, max_samples, 25):
            # create 100 data sets to test 
            for dataset_iter in xrange(0, 100):
                # create a dataset
                X1 = rng.randn(sample_size)
                X2 = rng.randn(sample_size)
                X3 = np.sqrt(1 - 2 * alpha_val**2) * rng.randn(sample_size) \
                    + alpha_val * (X1 + X2)
                X = np.c_[X1, X2, X3]
                y = np.dot(X, coef)
                
                # get the lasso's coefficients
                alphas, _, coefs = lars_path(X, y, method='lasso') 
                # get the non-negative garotte's coefficients
                ng_coefs, _, _ = non_negative_garotte_path(X, y, alpha_val)

                # test if either model's solution path matches the orinial model
                if np.any(np.all(ng_coefs.astype(np.bool) == coef.astype(np.bool))):
                    ng_path_correct = ng_path_correct + 1.0
                if np.any(np.all(coefs.astype(np.bool) == coef.astype(np.bool)[:, np.newaxis], axis=0)):
                    lars_path_correct = lars_path_correct + 1.0
            
            hits_pers_lars = lars_path_correct/100
            hits_lars.append(hits_pers_lars)
            lars_path_correct = 0  
            hits_pers_ng = ng_path_correct/100
            hits_ng.append(hits_pers_ng)
            ng_path_correct = 0
            
            
        pl.plot(xrange(25, max_samples, 25), hits_lars, 'r-')
        pl.plot(xrange(25, max_samples, 25), hits_ng, 'b-')
        pl.xlim([0, max_samples])
        pl.ylim([0, 1.1])
        hits_lars = []
        hits_ng = []

    pl.show()
