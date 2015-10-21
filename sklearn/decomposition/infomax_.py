"""
An Information Maximisation ICA (InfomaxICA) algorithm.

Reference: Based on the publications of Bell & Sejnowski 1995 (Infomax)
           and Lee, Girolami & Sejnowski, 1999 (extended Infomax)
"""
# Authors: Lukas Breuer <l.breuer@fz-juelich.de>
#          Juergen Dammers <j.dammers@fz-juelich.de>
#          Denis A. Engeman <denis.engemann@gemail.com>
#
# License: BSD (3-clause)

import math
import warnings
import numpy as np
from scipy import linalg
from scipy.stats import kurtosis

from ..base import BaseEstimator, TransformerMixin
from . import RandomizedPCA
from ..utils import check_random_state, check_array, as_float_array
from ..utils.extmath import fast_dot
from ..utils.validation import check_is_fitted
from ..utils.validation import FLOAT_DTYPES

__all__ = ['infomax', 'InfomaxICA']


def infomax(X, n_components=None, whiten=True, weights=None, l_rate=None,
            block=None, w_change=1e-12,
            anneal_deg=60., anneal_step=0.9, extended=False, n_subgauss=1,
            kurt_size=6000, ext_blocks=1, max_iter=200,
            random_state=None, blowup=1e4, blowup_fac=0.5, n_small_angle=20,
            use_bias=True, compute_sources=True, verbose=None):
    """Run the (extended) Infomax ICA decomposition on raw data

    Parameters
    ----------
    X : np.ndarray, shape (n_samples, n_features)
        The data to unmix.

    n_components : int, optional
        Number of components to extract. If None no dimension reduction
        is performed.

    whiten : boolean, optional
        If True perform an initial whitening of the data.
        If False, the data is assumed to have already been
        preprocessed: it should be centered, normed and white. Otherwise,
        you will get incorrect results.
        In this case, the parameter n_components will be ignored.
        Defaults to True.

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
        False.

    n_subgauss : int
        The number of subgaussian components. Only considered for extended
        Infomax.

    kurt_size : int
        The window size for kurtosis estimation. Only considered for extended
        Infomax.

    ext_blocks : int
        Only considered for extended Infomax.
        If positive, it denotes the number of blocks after which to recompute
        the Kurtosis, which is used to estimate the signs of the sources.
        In this case the number of sub-gaussian sources is automatically
        determined.
        If negative, the number of sub-gaussian sources to be used is fixed
        and equal to n_subgauss. In this case the Kurtosis is not estimated.

    max_iter : int
        The maximum number of iterations. Defaults to 200.

    random_state : int | np.random.RandomState
        If random_state is an int, use random_state as seed of the random
        number generator.
        If random_state is already a np.random.RandomState instance, use
        random_state as random number generator.

    blowup : float
        The maximum difference allowed between two succesive estimations of the
        unmixing matrix. Defaults to 1e4

    blowup_fac : float
        The factor by which the learning rate will be reduced if the
        difference between two succesive estimations of the
        unmixing matrix exceededs ``blowup``:
            l_rate *= blowup_fac
        Defaults to 0.5

    n_small_angle : int | None
        The maximum number of allowed steps in which the angle between two
        succesive estimations of the unmixing matrix is less than
        ``anneal_deg``.
        If None, this parameter is not taken into account to stop the
        iterations.
        Defaults to 20

    use_bias : bool
        This quantity indicates if the bias should be computed.
        Defaults to True

    compute_sources : boolean, optional
        If False, sources are not computed but only the rotation matrix.
        This can save memory when working with big data.
        Defaults to True.

    verbose : bool, str, int, or None
        If not None, override default verbose level.

    Returns
    -------
    unmixing_matrix : np.ndarray of float, shape (n_features, n_features)
        The linear unmixing operator.
    """
    rng = check_random_state(random_state)

    # define some default parameter
    max_weight = 1e8
    restart_fac = 0.9
    min_l_rate = 1e-10
    degconst = 180.0 / np.pi

    # for extended Infomax
    extmomentum = 0.5
    signsbias = 0.02
    signcount_threshold = 25
    signcount_step = 2

    # check data shape
    X = as_float_array(X, copy=False)  # to avoid casting problem
    n_samples, n_features = X.shape

    if not whiten and n_components is not None:
        n_components = None
        warnings.warn('Ignoring n_components with whiten=False.')

    if n_components is None:
        n_components = min(n_samples, n_features)
    if n_components > min(n_samples, n_features):
        n_components = min(n_samples, n_features)
        print("n_components is too large: it will be set to %s" % n_components)

    if whiten:
        pca = RandomizedPCA(n_components=n_components, whiten=True,
                            copy=True, random_state=random_state)
        X1 = pca.fit_transform(X)
        X_mean = pca.mean_
        prewhitening = pca.components_
        exp_var = pca.explained_variance_
        prewhitening *= np.sqrt(exp_var[:, None])
    else:
        # X must be casted to floats to avoid typing issues with numpy
        # 2.0 and the line below
        X1 = X  # copy has been taken care of
    n_features = n_components
    n_features_square = n_features ** 2

    # check input parameter
    # heuristic default - may need adjustment for
    # large or tiny data sets
    if l_rate is None:
        l_rate = 0.01 / math.log(n_features ** 2 + 1)

    if block is None:
        block = int(math.floor(math.sqrt(n_samples / 3.0)))

    if extended:
        print('computing%sInfomax ICA' % ' Extended ')

    if block == 0:
        block += 1
    # collect parameter
    nblock = n_samples // block
    lastt = (nblock - 1) * block + 1

    # initialize training
    if weights is None:
        # initialize weights as identity matrix
        weights = np.identity(n_features, dtype=np.float64)

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
    initial_ext_blocks = ext_blocks   # save the initial value in case of reset

    # for extended Infomax
    if extended:
        signs = np.ones(n_features)

        for k in range(n_subgauss):
            signs[k] = -1

        kurt_size = min(kurt_size, n_samples)
        old_kurt = np.zeros(n_features, dtype=np.float64)
        oldsigns = np.zeros(n_features)

    # trainings loop
    olddelta, oldchange = 1., 0.
    while step < max_iter:
        # shuffle data at each step
        permute = rng.permutation(n_samples)

        # ICA training block
        # loop across block samples
        for t in range(0, lastt, block):
            u = np.dot(X1[permute[t:t + block], :], weights)
            u += np.dot(bias, onesrow).T

            if extended:
                # extended ICA update
                y = np.tanh(u)
                weights += l_rate * np.dot(weights,
                                           BI -
                                           signs[None, :] * np.dot(u.T, y) -
                                           np.dot(u.T, u))
                if use_bias:
                    bias += l_rate * np.reshape(np.sum(y, axis=0,
                                                dtype=np.float64) * -2.0,
                                                (n_features, 1))

            else:
                # logistic ICA weights update
                y = 1.0 / (1.0 + np.exp(-u))
                weights += l_rate * np.dot(weights,
                                           BI + np.dot(u.T, (1.0 - 2.0 * y)))

                if use_bias:
                    bias += l_rate * np.reshape(np.sum((1.0 - 2.0 * y), axis=0,
                                                dtype=np.float64),
                                                (n_features, 1))

            # check change limit
            max_weight_val = np.max(np.abs(weights))
            if max_weight_val > max_weight:
                wts_blowup = True

            blockno += 1
            if wts_blowup:
                break

            # ICA kurtosis estimation
            if extended:

                if ext_blocks > 0 and blockno % ext_blocks == 0:

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
                    signs = np.sign(kurt + signsbias)

                    ndiff = (signs - oldsigns != 0).sum()
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
            if step > 2:
                angledelta = math.acos(np.sum(delta * olddelta) /
                                       math.sqrt(change * oldchange))
                angledelta *= degconst

            if verbose:
                print('step %d - lrate %5f, wchange %8.8f, angledelta %4.1f deg'
                      % (step, l_rate, change, angledelta))

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

                if n_small_angle is not None:
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

            ext_blocks = initial_ext_blocks

            # for extended Infomax
            if extended:
                signs = np.ones(n_features)
                for k in range(n_subgauss):
                    signs[k] = -1
                oldsigns = np.zeros(n_features)

            if l_rate > min_l_rate:
                if verbose:
                    print('... lowering learning rate to %g'
                          '\n... re-starting...' % l_rate)
            else:
                raise ValueError('Error in Infomax ICA: unmixing_matrix matrix'
                                 'might not be invertible!')

    del X1
    # prepare return values
    weights = weights.T
    if whiten is True:
        weights /= np.sqrt(exp_var)[None, :]
        if compute_sources:
            X -= X_mean
            sources = fast_dot(weights, fast_dot(prewhitening, X.T)).T
        else:
            sources = None
        return prewhitening, weights, sources, X_mean, step
    else:
        if compute_sources:
            sources = fast_dot(weights, X.T).T
        else:
            sources = None
        return None, weights, sources, None, step


class InfomaxICA(BaseEstimator, TransformerMixin):
    """An algorithm for Independent Component Analysis based on
    Information Maximisation.

    Parameters
    ----------
    n_components : int, optional
        Number of components to extract. If None is passed, all are used.

    whiten : boolean, optional
        If whiten is False, the data is already considered to be whitened,
        and no whitening is performed.

    w_init : np.ndarray, shape (n_features, n_features)
        The initialization of unmixing matrix. Defaults to None. If None, the
        identity matrix is used.

    l_rate : float
        This parameter indicates the relative size of the change in weights.
        Defaults to 0.01 / alog(n_features ^ 2.0)
        Note: Smaller learning rates will slow down the procedure.

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
        Whether to use the extended infomax algorithm or not. Defaults to False.

    n_subgauss : int
        The number of subgaussian components. Only considered if extended Infomax
        is True.

    kurt_size : int
        The window size for kurtosis estimation. Only considered if extended
        Infomax is True.

    ext_blocks : int
        Only considered if extended Infomax is True.
        If positive, it denotes the number of blocks after which to recompute
        the Kurtosis, which is used to estimate the signs of the sources.
        In this case the number of sub-gaussian sources is automatically
        determined.
        If negative, the number of sub-gaussian sources to be used is fixed
        and equal to n_subgauss. In this case, the Kutosis is not estimated.

    max_ter : int
        The maximum number of iterations. Defaults to 200.

    random_state : int | np.random.RandomState
        If random_state is an int, use random_state as seed of the rng.
        If random_state is already a np.random.RandomState instance, use
        random_state as rng.

    blowup : float
        The maximum difference allowed between two successive estimations of the
        unmixing matrix. Defaults to 1e4

    blowup_fac : float
        The factor by which the learning rate will be reduced if the difference
        between two successive estimations of the unmixing matrix exceeds
        ``blowup``:
            l_rate *= blowup_fac
        Defaults to 0.5

    n_small_angle : int | None
        The maximum number of allowed steps in which the angle between two
        successive estimations of the unmixing matrix is less than ``anneal_deg``.
        If None, this parameter is not taken into account to stop the iterations.
        Defaults to 20

    use_bias : boolean, optional
        This parameter indicates if the bias should be computed.
        Defaults to True.

    Attributes
    ----------
    """
    def __init__(self, n_components=None, whiten=True, weights=None, l_rate=None,
                 block=None, w_change=1e-12, anneal_deg=60., anneal_step=0.9,
                 extended=False, n_subgauss=1, kurt_size=6000, ext_blocks=1,
                 max_iter=200, random_state=None, blowup=1e4, blowup_fac=0.5,
                 n_small_angle=20, use_bias=True, verbose=None):
        super(InfomaxICA, self).__init__()
        self.n_components = n_components
        self.whiten = whiten
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
        self.blowup = blowup
        self.blowup_fac = blowup_fac
        self.n_small_angle = n_small_angle
        self.use_bias = use_bias

    def _fit(self, X, compute_sources=False):
        """Fit the model

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        compute_sources : boolean
            If False, sources are not computed but only rotation matrix.
            This can save memory when working with big data.
            Defaults to False.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)

        """
        whitening, unmixing, sources, X_mean, n_iter = infomax(
            X=X, n_components=self.n_components, whiten=self.whiten,
            weights=self.weights, l_rate=self.l_rate, block=self.block,
            w_change=self.w_change, anneal_deg=self.anneal_deg,
            anneal_step=self.anneal_step, extended=self.extended,
            n_subgauss=self.n_subgauss, kurt_size=self.kurt_size,
            ext_blocks=self.ext_blocks, max_iter=self.max_iter,
            random_state=self.random_state, blowup=self.blowup,
            blowup_fac=self.blowup_fac, n_small_angle=self.n_small_angle,
            use_bias=self.use_bias, compute_sources=compute_sources)

        if self.whiten:
            self.components_ = np.dot(unmixing, whitening)
            self.mean_ = X_mean
            self.whitening_ = whitening
        else:
            self.components_ = unmixing

        self.mixing_ = linalg.pinv(self.components_)
        self.n_iter_ = n_iter

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
        self._fit(X, compute_sources=True)
        return self

    def transform(self, X, y=None, copy=True):
        """Recover the sources from X (apply the unmixing matrix).

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Data to transform, where n_samples is the number of samples
            and n_features are the number of features.

        copy : boolean, optional
            If False, data passed to fit are overwritten.
            Defaults to True.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)

        """
        check_is_fitted(self, 'mixing_')

        X = check_array(X, copy=copy, dtype=FLOAT_DTYPES)
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
        copy : boolean, optional
            If False, data passed to fit are overwritten.
            Defaults to True.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_features)
        """
        check_is_fitted(self, 'mixing_')

        X = check_array(X, copy=(copy and self.whiten), dtype=FLOAT_DTYPES)
        X = fast_dot(X, self.mixing_.T)
        if self.whiten:
            X += self.mean_

        return X
