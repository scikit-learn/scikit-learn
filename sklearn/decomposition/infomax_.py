# Authors: Lukas Breuer <l.breuer@fz-juelich.de>
#          Juergen Dammers <j.dammers@fz-juelich.de>
#          Denis A. Engeman <denis.engemann@gemail.com>
#
# License: BSD (3-clause)

import warnings
import numpy as np
from scipy import linalg
import math
from scipy.stats import kurtosis

from ..base import BaseEstimator, TransformerMixin
from ..externals import six
from ..externals.six import moves
from ..utils import array2d, as_float_array, check_random_state
from ..utils.extmath import fast_dot


def infomax(X, n_components=None, weights=None, l_rate=None, block=None,
            w_change=1e-12, anneal_deg=60., anneal_step=0.9, extended=False,
            n_subgauss=1, kurt_size=6000, ext_blocks=1, max_iter=200,
            random_state=None, whiten=True, return_X_mean=False,
            compute_sources=True):
    """Run the (extended) Infomax ICA decomposition on raw data

    based on the publications of Bell & Sejnowski 1995 (Infomax)
    and Lee, Girolami & Sejnowski, 1999 (extended Infomax)

    Parameters
    ----------
    X : np.ndarray, shape (n_samples, n_features)
        The data to unmix.

    n_components : int, optional
        Number of components to extract. If None no dimension reduction
        is performed.

    w_init : np.ndarray, shape (n_features, n_features)
        The initialized unmixing matrix. Defaults to None. If None, the
        identity matrix is used.

    l_rate : float
        This quantity indicates the relative size of the change in weights.
        Note. Smaller learining rates will slow down the procedure.
        Defaults to 0.010d / alog(n_features ^ 2.0)

    block : int
        The block size of randomly chosen data segment.
        Defaults to floor(sqrt(n_times / 3d))

    w_change : float
        The change at which to stop iteration. Defaults to 1e-12.

    anneal_deg : float
        The angle at which (in degree) the learning rate will be reduced.
        Defaults to 60.0

    anneal_step : float
        The factor by which the learning rate will be reduced once
        ``anneal_deg`` is exceeded:
            l_rate *= anneal_step
        Defaults to 0.9

    extended : bool
        Wheather to use the extended infomax algorithm or not. Defaults to
        True.

    n_subgauss : int
        The number of subgaussian components. Only considered for extended
        Infomax.

    kurt_size : int
        The window size for kurtosis estimation. Only considered for extended
        Infomax.

    ext_blocks : int
        The number of blocks after which to recompute Kurtosis.
        Only considered for extended Infomax.

    max_iter : int
        The maximum number of iterations. Defaults to 200.

    whiten : boolean, optional
        If True perform an initial whitening of the data.
        If False, the data is assumed to have already been
        preprocessed: it should be centered, normed and white.
        Otherwise you will get incorrect results.
        In this case the parameter n_components will be ignored.

    return_X_mean : bool, optional
        If True, X_mean is returned too.

    compute_sources : bool, optional
        If False, sources are not computed, but only the rotation matrix.
        This can save memory when working with big data. Defaults to True.


    Returns
    -------
    prewhitening : array, shape (n_components, n_features) | None.
        If whiten is 'True', K is the pre-whitening matrix that projects data
        onto the first n_components principal components. If whiten is 'False',
        K is 'None'.

    weights : array, shape (n_components, n_components)
        Estimated un-mixing matrix.
        The mixing matrix can be obtained by::

            w = np.dot(W, K.T)
            A = w.T * (w * w.T).I

    sources : array, shape (n_components, n_samples) | None
        Estimated source matrix

    X_mean : array, shape (n_features, )
        The mean over features. Returned only if return_X_mean is True.

    Notes
    -----

    The data matrix X is considered to be a linear combination of
    non-Gaussian (independent) components i.e. X = AS where columns of S
    contain the independent components and A is a linear mixing
    matrix. In short ICA attempts to `un-mix' the data by estimating an
    un-mixing matrix W where `S = W K X.``
    """
    rng = check_random_state(random_state)

    # define some default parameters (constant-like, not exposed for now)
    max_weight = 1e8
    restart_fac = 0.9
    min_l_rate = 1e-10
    blowup = 1e4
    blowup_fac = 0.5
    n_small_angle = 20
    degconst = 180.0 / np.pi

    # for extended Infomax (constant-like, not exposed for now)
    extmomentum = 0.5
    signsbias = 0.02
    signcount_threshold = 25
    signcount_step = 2
    if ext_blocks > 0:  # allow not to recompute kurtosis
        n_subgauss = 1  # but initialize n_subgauss to 1 if you recompute

    # check data shape
    n_samples, n_features = X.shape
    X = array2d(X, copy=whiten).T

    if not whiten and n_components is not None:
        n_components = None
        warnings.warn('Ignoring n_components with whiten=False.')

    if n_components is None:
        n_components = min(n_samples, n_features)
    if n_components > min(n_samples, n_features):
        n_components = min(n_samples, n_features)
        print("n_components is too large: it will be set to %s" % n_components)

    if whiten:
        # Centering the columns (ie the variables)
        X_mean = X.mean(axis=-1)
        X -= X_mean[:, np.newaxis]

        # Whitening and preprocessing by PCA
        u, d, _ = linalg.svd(X, full_matrices=False)

        del _
        prewhitening = (u / d).T[:n_components]  # see (6.33) p.140
        del u, d
        X1 = np.dot(prewhitening, X)
        # see (13.6) p.267 Here X1 is white and data
        # in X has been projected onto a subspace by PCA
        X1 *= np.sqrt(n_features)
    else:
        # X must be casted to floats to avoid typing issues with numpy
        # 2.0 and the line below
        X1 = as_float_array(X, copy=False)  # copy has been taken care of
    X1 = X1.T
    n_features = n_components
    n_features_square = n_features ** 2
    # check input parameter
    # heuristic default - may need adjustment for
    # large or tiny X sets
    if l_rate is None:
        l_rate = 0.01 / math.log(n_features ** 2.0)

    if block is None:
        block = int(math.floor(math.sqrt(n_samples / 3.0)))

    # collect parameter
    nblock = n_samples // block
    lastt = (nblock - 1) * block + 1

    # initialize training
    if weights is None:
        # initialize weights as identity matrix
        weights = np.identity(n_components, dtype=np.float64)

    BI = block * np.identity(n_features, dtype=np.float64)
    bias = np.zeros((n_features, 1), dtype=np.float64)
    onesrow = np.ones((1, block), dtype=np.float64)
    startweights = weights.copy()
    oldweights = startweights.copy()
    step = 0
    count_small_angle = 0
    wts_blowup = False
    blockno = 0
    signcount = 0

    # for extended Infomax
    if extended is True:
        signs = np.identity(n_features)
        signs.flat[slice(0, n_features * n_subgauss, n_features)]
        kurt_size = min(kurt_size, n_samples)
        old_kurt = np.zeros(n_features, dtype=np.float64)
        oldsigns = np.zeros((n_features, n_features))

    # trainings loop
    olddelta, oldchange = 1., 0.
    while step < max_iter:

        # shuffle X at each step
        rng.seed(step)  # --> permutation is fixed but differs at each step
        permute = list(range(n_samples))
        rng.shuffle(permute)

        # ICA training block
        # loop across block samples
        for t in range(0, lastt, block):
            u = np.dot(X1[permute[t:t + block], :], weights)
            u += np.dot(bias, onesrow).T

            if extended is True:
                # extended ICA update
                y = np.tanh(u)
                weights += l_rate * np.dot(weights,
                                           BI - np.dot(np.dot(u.T, y), signs) -
                                           np.dot(u.T, u))
                bias += l_rate * np.reshape(np.sum(y, axis=0,
                                            dtype=np.float64) * -2.0,
                                            (n_features, 1))

            else:
                # logistic ICA weights update
                y = 1.0 / (1.0 + np.exp(-u))
                weights += l_rate * np.dot(weights,
                                           BI + np.dot(u.T, (1.0 - 2.0 * y)))
                bias += l_rate * np.reshape(np.sum((1.0 - 2.0 * y), axis=0,
                                            dtype=np.float64), (n_features, 1))

            # check change limit
            max_weight_val = np.max(np.abs(weights))
            if max_weight_val > max_weight:
                wts_blowup = True

            blockno += 1
            if wts_blowup:
                break

            # ICA kurtosis estimation
            if extended is True:

                n = np.fix(blockno / ext_blocks)

                if np.abs(n) * ext_blocks == blockno:
                    if kurt_size < n_samples:
                        rp = np.floor(rng.uniform(0, 1, kurt_size) *
                                      (n_samples - 1))
                        tpartact = np.dot(X1[rp.astype(int), :], weights).T
                    else:
                        tpartact = np.dot(X1, weights).T

                    # estimate kurtosis
                    kurt = kurtosis(tpartact, axis=1, fisher=True)

                    if extmomentum != 0:
                        kurt = (extmomentum * old_kurt +
                                (1.0 - extmomentum) * kurt)
                        old_kurt = kurt

                    # estimate weighted signs
                    signs.flat[::n_features + 1] = ((kurt + signsbias) /
                                                    np.abs(kurt + signsbias))

                    ndiff = ((signs.flat[::n_features + 1] -
                              oldsigns.flat[::n_features + 1]) != 0).sum()
                    if ndiff == 0:
                        signcount += 1
                    else:
                        signcount = 0
                    oldsigns = signs

                    if signcount >= signcount_threshold:
                        ext_blocks = np.fix(ext_blocks * signcount_step)
                        signcount = 0

        # here we continue after the for
        # loop over the ICA training blocks
        # if weights in bounds:
        if not wts_blowup:
            oldwtchange = weights - oldweights
            step += 1
            angledelta = 0.0
            delta = oldwtchange.reshape(1, n_features_square)
            change = np.sum(delta * delta, dtype=np.float64)
            if step > 1:
                angledelta = math.acos(np.sum(delta * olddelta) /
                                       math.sqrt(change * oldchange))
                angledelta *= degconst

            # anneal learning rate
            oldweights = weights.copy()
            if angledelta > anneal_deg:
                l_rate *= anneal_step    # anneal learning rate
                # accumulate angledelta until anneal_deg reached l_rates
                olddelta = delta
                oldchange = change
                count_small_angle = 0  # reset count when angle delta is large
            else:
                if step == 1:  # on first step only
                    olddelta = delta  # initialize
                    oldchange = change
                count_small_angle += 1
                if count_small_angle > n_small_angle:
                    max_iter = step

            # apply stopping rule
            if step > 2 and change < w_change:
                step = max_iter
            elif change > blowup:
                l_rate *= blowup_fac

        # restart if weights blow up
        # (for lowering l_rate)
        else:
            step = 0  # start again
            wts_blowup = 0  # re-initialize variables
            blockno = 1
            l_rate *= restart_fac  # with lower learning rate
            weights = startweights.copy()
            oldweights = startweights.copy()
            olddelta = np.zeros((1, n_features_square), dtype=np.float64)
            bias = np.zeros((n_features, 1), dtype=np.float64)

            # for extended Infomax
            if extended:
                signs = np.identity(n_features)
                signs.flat[slice(0, n_features * n_subgauss, n_features)]
                oldsigns = np.zeros((n_features, n_features))

            if l_rate > min_l_rate:
                pass  # maybe add some logging later
            else:
                raise ValueError('Error in Infomax ICA: unmixing_matrix matrix'
                                 'might not be invertible!')

    del X1
    weights = weights.T
    if whiten:
        if compute_sources:
            sources = fast_dot(fast_dot(weights, prewhitening), X).T
        else:
            sources = None
        if return_X_mean:
            return prewhitening, weights, sources, X_mean
        else:
            return prewhitening, weights, sources
    else:
        if compute_sources:
            sources = fast_dot(weights, X).T
        else:
            sources = None
        if return_X_mean:
            return None, weights, sources, None
        else:
            return None, weights, sources


class InfomaxICA(BaseEstimator, TransformerMixin):
    """Run the (extended) Infomax ICA decomposition on raw data

    based on the publications of Bell & Sejnowski 1995 (Infomax)
    and Lee, Girolami & Sejnowski, 1999 (extended Infomax)

    Parameters
    ----------
    n_components : int, optional
        Number of components to extract. If None no dimension reduction
        is performed.

    w_init : np.ndarray, shape (n_features, n_features)
        The initialized unmixing matrix. Defaults to None. If None, the
        identity matrix is used.

    l_rate : float
        This quantity indicates the relative size of the change in weights.
        Note. Smaller learining rates will slow down the procedure.
        Defaults to 0.010d / alog(n_features ^ 2.0)

    block : int
        The block size of randomly chosen data segment.
        Defaults to floor(sqrt(n_times / 3d))

    w_change : float
        The change at which to stop iteration. Defaults to 1e-12.

    anneal_deg : float
        The angle at which (in degree) the learning rate will be reduced.
        Defaults to 60.0

    anneal_step : float
        The factor by which the learning rate will be reduced once
        ``anneal_deg`` is exceeded:
            l_rate *= anneal_step
        Defaults to 0.9

    extended : bool
        Wheather to use the extended infomax algorithm or not. Defaults to
        True.

    n_subgauss : int
        The number of subgaussian components. Only considered for extended
        Infomax.

    kurt_size : int
        The window size for kurtosis estimation. Only considered for extended
        Infomax.

    ext_blocks : int
        The number of blocks after which to recompute Kurtosis.
        Only considered for extended Infomax.

    max_iter : int
        The maximum number of iterations. Defaults to 200.

    whiten : boolean, optional
        If True perform an initial whitening of the data.
        If False, the data is assumed to have already been
        preprocessed: it should be centered, normed and white.
        Otherwise you will get incorrect results.
        In this case the parameter n_components will be ignored.
    """
    def __init__(self, n_components=None, weights=None, l_rate=None,
                 block=None, w_change=1e-12, anneal_deg=60., anneal_step=0.9,
                 extended=False, n_subgauss=1, kurt_size=6000, ext_blocks=1,
                 max_iter=200, random_state=None, whiten=True):
        self.n_components = n_components
        self.weights = weights
        self.l_rate = l_rate
        self.block = block
        self.w_change = w_change
        self.anneal_deg = anneal_deg
        self.anneal_step = anneal_step
        self.extended = extended
        self.n_subgauss = n_subgauss
        self.kurt_size = kurt_size
        self.ext_blocks = ext_blocks
        self.max_iter = max_iter
        self.random_state = random_state
        self.whiten = whiten

    def _fit(self, X, compute_sources=False):
        """Fit the model

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        compute_sources : bool
            If False, sources are not computes but only the rotation matrix.
            This can save memory when working with big data. Defaults to False.

        Returns
        -------
            X_new : array-like, shape (n_samples, n_components)
        """
        whitening, unmixing, sources, X_mean = infomax(
            X=X, n_components=self.n_components, weights=self.weights,
            l_rate=self.l_rate, block=self.block, w_change=self.w_change,
            anneal_deg=self.anneal_deg, anneal_step=self.anneal_step,
            extended=self.extended, n_subgauss=self.n_subgauss,
            kurt_size=self.kurt_size, ext_blocks=self.ext_blocks,
            max_iter=self.max_iter, random_state=self.random_state,
            whiten=self.whiten, return_X_mean=True,
            compute_sources=compute_sources)

        if self.whiten:
            self.components_ = np.dot(unmixing, whitening)
            self.mean_ = X_mean
            self.whitening_ = whitening
        else:
            self.components_ = unmixing

        self.mixing_ = linalg.pinv(self.components_)

        if compute_sources:
            self.__sources = sources

        return sources

    def fit_transform(self, X, y=None):
        """Fit the model and recover the sources from X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)
        """
        return self._fit(X, compute_sources=True)

    def fit(self, X, y=None):
        """Fit the model to X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        self
        """
        self._fit(X, compute_sources=True)  # will become False in 0.16
        return self

    def transform(self, X, y=None, copy=True):
        """Recover the sources from X (apply the unmixing matrix).

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Data to transform, where n_samples is the number of samples
            and n_features is the number of features.

        copy : bool (optional)
            If False, data passed to fit are overwritten. Defaults to True.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)
        """
        X = array2d(X, copy=copy)
        if self.whiten:
            X -= self.mean_

        return fast_dot(X, self.components_.T)

    def inverse_transform(self, X, copy=True):
        """Transform the sources back to the mixed data (apply mixing matrix).

        Parameters
        ----------
        X : array-like, shape (n_samples, n_components)
            Sources, where n_samples is the number of samples
            and n_components is the number of components.
        copy : bool (optional)
            If False, data passed to fit are overwritten. Defaults to True.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_features)
        """
        if copy:
            X = X.copy()
        X = fast_dot(X, self.mixing_.T)
        if self.whiten:
            X += self.mean_

        return X
