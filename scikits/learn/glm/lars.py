
import numpy as np
from scipy import linalg
from .base import LinearModel
import scipy.sparse as sp # needed by LeastAngleRegression
from .._minilearn import lars_fit_wrap


# Notes: np.ma.dot copies the masked array before doing the dot
# product. Maybe we should implement in C our own masked_dot that does
# not make unnecessary copies.

# all linalg.solve solve a triangular system, so this could be heavily
# optimized by binding (in scipy ?) trsv or trsm

def lars(X, y, max_iter=None, alpha_min=0, method="lar", precompute=True):
    """
    lar -> m <= x.shape[1]
    lasso -> m can be > x.shape[1]

    precompute : empty for now

    TODO: detect stationary points.
    Lasso variant
    store full path
    """

    n_samples, n_features = X.shape

    if max_iter is None:
        max_iter = min(n_samples, n_features)

    max_pred = max_iter # OK for now

    mu       = np.zeros (X.shape[0])
    beta     = np.zeros ((max_iter + 1, X.shape[1]))
    alphas   = np.zeros (max_iter + 1)
    n_iter, n_pred = 0, 0
    active   = list()
    # holds the sign of covariance
    sign_active = np.empty (max_pred, dtype=np.int8)
    drop = False

    # will hold the cholesky factorization
    # only lower part is referenced. We do not create it as
    # empty array because chol_solve calls chkfinite on the
    # whole array, which can cause problems.
    L = np.zeros ((max_pred, max_pred), dtype=np.float64)

    Xt  = X.T
    Xna = Xt.view(np.ma.MaskedArray) # variables not in the active set
                                     # should have a better name

    while 1:


        # Calculate covariance matrix and get maximum
        res = y - np.dot (X, beta[n_iter]) # there are better ways
        Cov = np.ma.dot (Xna, res)

        imax    = np.ma.argmax (np.ma.abs(Cov), fill_value=0.) #rename
        Cov_max =  (Cov [imax])
        Cov[imax] = np.ma.masked


        alpha = np.abs(Cov_max) #sum (np.abs(beta[n_iter]))
        alphas [n_iter] = np.max(np.abs(np.dot(Xt, res))) #sum (np.abs(beta[n_iter]))
        if (n_iter >= max_iter or n_pred >= max_pred ):
            break

        if (alpha < alpha_min): break


        if not drop:

            # Update the Cholesky factorization of (Xa * Xa') #
            #                                                 #
            #          ( L   0 )                              #
            #   L  ->  (       )  , where L * w = b           #
            #          ( w   z )    z = 1 - ||w||             #
            #                                                 #
            #   where u is the last added to the active set   #

            n_pred += 1
            active.append(imax)
            Xna[imax] = np.ma.masked
            sign_active [n_pred-1] = np.sign (Cov_max)

            X_max = Xt[imax]

            c = np.dot (X_max, X_max)
            L [n_pred-1, n_pred-1] = c

            if n_pred > 1:
                b = np.dot (X_max, Xa.T)

                # please refactor me, using linalg.solve is overkill
                L [n_pred-1, :n_pred-1] = linalg.solve (L[:n_pred-1, :n_pred-1], b)
                v = np.dot(L [n_pred-1, :n_pred-1], L [n_pred - 1, :n_pred -1])
                L [n_pred-1,  n_pred-1] = np.sqrt (c - v)
        else:
            drop = False

        Xa = Xt[active] # also Xna[~Xna.mask]

        # Now we go into the normal equations dance.
        # (Golub & Van Loan, 1996)

        b = np.copysign (Cov_max.repeat(n_pred), sign_active[:n_pred])
        b = linalg.cho_solve ((L[:n_pred, :n_pred], True),  b)

        C = A = np.abs(Cov_max)
        u = np.dot (Xa.T, b)
        a = np.ma.dot (Xna, u)

        # equation 2.13, there's probably a simpler way
        g1 = (C - Cov) / (A - a)
        g2 = (C + Cov) / (A + a)

        # one for the border cases
        g = np.ma.concatenate((g1, g2))

        g = g[g >= 0.]
        gamma_ = np.ma.min (g)

        if method == 'lasso':

            z = - beta[n_iter, active] / b
            z[z <= 0.] = np.inf

            idx = np.argmin(z)

            if z[idx] < gamma_:
                gamma_ = z[idx]
                drop = True

        n_iter += 1
        beta[n_iter, active] = beta[n_iter - 1, active] + gamma_ * b
        if n_iter == X.shape[1]:
            # oveeeerkilll 
            beta[n_iter] =  linalg.lstsq (X, y)[0]
            alpha = 0.
            break

        if drop:
            n_pred -= 1
            drop_idx = active.pop (idx)
            print 'dropped ', idx, ' at ', n_iter, ' iteration'
            Xa = Xt[active] # duplicate
            L[:n_pred, :n_pred] = linalg.cholesky(np.dot(Xa, Xa.T), lower=True)
            sign_active = np.delete(sign_active, idx) # do an append to maintain size
            # should be done using cholesky deletes

    if alpha < alpha_min: # interpolate
        # interpolation factor 0 <= ss < 1
        ss = (alphas[n_iter-1] - alpha_min) / (alphas[n_iter-1] - alphas[n_iter])
        beta[n_iter] = beta[n_iter-1] + ss*(beta[n_iter] - beta[n_iter-1]);
        alphas[n_iter] = alpha_min
        alphas = alphas[:n_iter+1]
        beta = beta[:n_iter+1]

    return alphas, active, beta.T


class LARS (LinearModel):

    def __init__(self, n_features, normalize=True):
        self.n_features = n_features
        self.normalize = normalize
        self.coef_ = None

    def fit (self, X, Y):
                # will only normalize non-zero columns

        if self.normalize:
            self._xmean = X.mean(0)
            self._ymean = Y.mean(0)
            X = X - self._xmean
            Y = Y - self._ymean
            self._norms = np.apply_along_axis (np.linalg.norm, 0, X)
            nonzeros = np.flatnonzero(self._norms)
            X[:, nonzeros] /= self._norms[nonzeros]

        method = 'lar'
        alphas_, active, coef_path_ = lars (X, Y,
                                max_iter=self.n_features, method=method)
        print alphas_
        self.coef_ = coef_path_[:,-1]
        return self


class LassoLARS (LinearModel):

    def __init__(self, alpha, normalize=True):
        self.alpha = alpha
        self.normalize = normalize
        self.coef_ = None

    def fit (self, X, Y):
                # will only normalize non-zero columns

        n_samples = X.shape[0]
        alpha = self.alpha * n_samples # scale alpha with number of samples

        if self.normalize:
            self._xmean = X.mean(0)
            self._ymean = Y.mean(0)
            X = X - self._xmean
            Y = Y - self._ymean
            self._norms = np.apply_along_axis (np.linalg.norm, 0, X)
            nonzeros = np.flatnonzero(self._norms)
            X[:, nonzeros] /= self._norms[nonzeros]

        method = 'lasso'
        alphas_, active, coef_path_ = lars (X, Y,
                                            alpha_min=alpha, method=method)
        self.coef_ = coef_path_[:,-1]
        return self


#### OLD C-based LARS : will probably be removed


class LeastAngleRegression(LinearModel):
    """
    Least Angle Regression using the LARS algorithm.

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        parameter vector (w in the fomulation formula)

    `intercept_` : float
        independent term in decision function.

    `coef_path_` : array, shape = [max_features + 1, n_features]
         Full coeffients path.

    Notes
    -----
    predict does only work correctly in the case of normalized
    predictors.

    See also
    --------
    scikits.learn.glm.Lasso

    """

    def __init__(self):
        self.alphas_ = np.empty(0, dtype=np.float64)
        self._chol   = np.empty(0, dtype=np.float64)
        self.beta_    = np.empty(0, dtype=np.float64)

    def fit (self, X, Y, fit_intercept=True, max_features=None, normalize=True):
        """
        Fit the model according to data X, Y.

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]
            Training data

        Y : numpy array of shape [n_samples]
            Target values

        fit_intercept : boolean, optional
            wether to calculate the intercept for this model. If set
            to false, no intercept will be used in calculations
            (e.g. data is expected to be already centered).

        max_features : int, optional
            number of features to get into the model. The iterative
            will stop just before the `max_features` variable enters
            in the active set. If not specified, min(N, p) - 1
            will be used.

        normalize : boolean
            whether to normalize (make all non-zero columns have mean
            0 and norm 1).
        """
        ## TODO: resize (not create) arrays, check shape,
        ##    add a real intercept

        X  = np.asanyarray(X, dtype=np.float64, order='C')
        _Y = np.asanyarray(Y, dtype=np.float64, order='C')

        if Y is _Y: Y = _Y.copy()
        else: Y = _Y

        if max_features is None:
            max_features = min(*X.shape)-1

        sum_k = max_features * (max_features + 1) /2
        self.alphas_.resize(max_features + 1)
        self._chol.resize(sum_k)
        self.beta_.resize(sum_k)
        coef_row = np.zeros(sum_k, dtype=np.int32)
        coef_col = np.zeros(sum_k, dtype=np.int32)


        if normalize:
            # will only normalize non-zero columns
            self._xmean = X.mean(0)
            self._ymean = Y.mean(0)
            X = X - self._xmean
            Y = Y - self._ymean
            self._norms = np.apply_along_axis (np.linalg.norm, 0, X)
            nonzeros = np.flatnonzero(self._norms)
            X[:, nonzeros] /= self._norms[nonzeros]
        else:
            self._xmean = 0.
            self._ymean = 0.

        lars_fit_wrap(0, X, Y, self.beta_, self.alphas_, coef_row,
                      coef_col, self._chol, max_features)

        self.coef_path_ = sp.coo_matrix((self.beta_,
                                        (coef_row, coef_col)),
                                        shape=(X.shape[1], max_features+1)).todense()

        self.coef_ = np.ravel(self.coef_path_[:, max_features])

        if fit_intercept:
            self.intercept_ = self._ymean
        else:
            self.intercept_ = 0.

        return self


    def predict(self, X, normalize=True):
        """
        Predict using the linear model.

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]

        Returns
        -------
        C : array, shape = [n_samples]
            Returns predicted values.
        """
        X = np.asanyarray(X, dtype=np.float64, order='C')
        if normalize:
            X -= self._xmean
            X /= self._norms
        return  np.dot(X, self.coef_) + self.intercept_


    
