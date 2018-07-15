
"""
Python implementation of the ICA algorithms: FastICA including, infomax, extendented infomax and picard.

    References
    ----------
    .. [1] Hyvärinen, A., 1999. Fast and robust fixed-point algorithms for
           independent component analysis. IEEE transactions on Neural
           Networks, 10(3), pp.626-634.
    .. [2] Bell, A.J., Sejnowski, T.J., 1995. An information-maximization
           approach to blind separation and blind deconvolution. Neural
           computation, 7(6), pp.1129-1159.
    .. [3] Lee, T.W., Girolami, M., Sejnowski, T.J., 1999. Independent
           component analysis using an extended infomax algorithm for mixed
           subgaussian and supergaussian sources. Neural computation, 11(2),
           pp.417-441.
    .. [4] Ablin, P., Cardoso, J.F., Gramfort, A., 2017. Faster Independent
           Component Analysis by preconditioning with Hessian approximations.
           arXiv:1706.08171
"""
# Authors: Nikesh Bajaj, with help of reference articles and mne library


import math
import numpy as np
from copy import deepcopy
from distutils.version import LooseVersion

class ICA():
    '''
       S =  A*X
       X =  W*S
       X : input data shape (nf,ns), nf- number of features or number of channels, ns- number of samples
       S : decomposed data  (n,ns)   n - number of components choosen, default n=nf
       A : Transform matirx  (n,n)
       W : inverse transform matrix (n,n)
    
    Signal decomposition using Independent Component Analysis (ICA), very usefule for EEG signal decompositions
    Including InfoMax, Extendent InfoMax and Picard methods, default as FastICA as usual 

    Parameters
    ----------
    n_components : int | None
        The number of components used for ICA decomposition. it must be smaller than 'max_pca_components'. 
        If None, all PCA components will be used default None, set to max number of components
    max_pca_components : int | None 
        The number of components used for PCA decomposition. If None, no dimensionality reduction will be 
        applied and `max_pca_components` will equal the number of channels (number of features) supplied for 
        decomposing data.
    n_pca_components:  int | float
        The number of PCA components used after ICA recomposition.
    random_state:  None | int | instance of np.random.RandomState
    method : {'fastica', 'infomax', 'extended-infomax', 'picard'}
        The ICA method to use. Defaults to 'fastica'. For reference, see [1]_,
        [2]_, [3]_ and [4]_. 
    fit_params : dict | None
        Additional parameters passed to the ICA estimator as specified by `method`.
    max_iter : int
        Maximum number of iterations during fit.
    
    
    Attributes
    ----------
    Estimated Values
    -----
    pca_mean_        :  mean substacted from data before computing PCA
    pca_components_  : PCA transform matrix
    pca_explained_variance_ :  variance of Principle components
    unmixing_matrix_ : ICA unmixing matrix A
    mixing_matrix_   : ICA mixing matrix W
    whitener_        : Standard deviaation of data before applying ICA 
    
    n_components
    max_pca_components
    n_pca_components
    random_state
    fit_params
    
    
    Methods:
    ---------
       S =  A*X
       X =  W*S
       X : input data shape (nf,ns), nf- number of features or number of channels, ns- number of samples
       S : decomposed data  (n,ns)   n - number of components choosen, default n=nf
       A : Transform matirx  (n,n)
       W : inverse transform matrix (n,n)
       
    fit(self, X, normalize=False):
        Fitting to data matrix X, X ndarray (nf,ns)
    transform(self, Xdata):
        Decompose Xdata into Independent Components
        return Xd (ndarray)
    get_tMatrix(self):
        Get Tranformation matrix
        return A (n,n)
    get_sMatrix(self):
        Get Inverse Transform matrix
        return W (n,n)
    whitening(self, X):
        To normlize the standard deviation of entire data (not the usual normailization)
    
    '''
    
       
    
    
    def __init__(self, n_components=None, max_pca_components=None, n_pca_components=None, random_state=None, 
        method='fastica', fit_params=None, max_iter=200):



        methods = ('fastica', 'infomax', 'extended-infomax', 'picard')
        if method not in methods:
            raise ValueError('method must be "%s". You passed: "%s"' % ('" or "'.join(methods), method))

        if (n_components is not None and max_pca_components is not None and n_components > max_pca_components):
            raise ValueError('n_components must be smaller than max_pca_components')
            
        self.n_components = n_components
        self.max_pca_components = max_pca_components   
        self.n_pca_components = n_pca_components
        self.random_state = random_state

        if fit_params is None:
            fit_params = {}
        fit_params = deepcopy(fit_params)  # avoid side effects
        if "extended" in fit_params:
            raise ValueError("'extended' parameter provided. You should "
                             "rather use method='extended-infomax'.")
        if method == 'fastica':
            update = {'algorithm': 'parallel', 'fun': 'logcosh',
                      'fun_args': None}
            fit_params.update(dict((k, v) for k, v in update.items() if k
                              not in fit_params))
        elif method == 'infomax':
            fit_params.update({'extended': False})
        elif method == 'extended-infomax':
            fit_params.update({'extended': True})
        elif method == 'picard':
            update = {'ortho': True, 'fun': 'tanh', 'tol': 1e-5}
            fit_params.update(dict((k, v) for k, v in update.items() if k
                              not in fit_params))
        if 'max_iter' not in fit_params:
            fit_params['max_iter'] = max_iter
        self.max_iter = max_iter
        self.fit_params = fit_params


        self.method = method
        
    def fit(self, X, normalize=False):
        """Run the ICA decomposition on X.
        
        
        X = array like: Shape (nf,ns) or (nCh, nSamples)
        
        """
        
        
        if self.max_pca_components is None:
            self.max_pca_components = X.shape[0]

        self.n_samples_ = X.shape[1]


        Xw, self.whitener_ = self.whitening(X)



        from sklearn.decomposition import PCA

        if not check_version('sklearn', '0.18'):
            pca = PCA(n_components=self.max_pca_components, whiten=True, copy=True)
        else:
            pca = PCA(n_components=self.max_pca_components, whiten=True, copy=True,
                      svd_solver='full')

        Xpca = pca.fit_transform(Xw.T)


        self.pca_mean_ = pca.mean_
        self.pca_components_ = pca.components_
        self.pca_explained_variance_ = exp_var = pca.explained_variance_
        if not check_version('sklearn', '0.16'):
            # sklearn < 0.16 did not apply whitening to the components, so we
            # need to do this manually
            self.pca_components_ *= np.sqrt(exp_var[:, None])
        del pca


        if self.method == 'fastica':
        	from sklearn.decomposition import FastICA
        	ica = FastICA(whiten=False, random_state=self.random_state, **self.fit_params)
        	ica.fit(Xpca)
        	self.unmixing_matrix_ = ica.components_

        elif self.method in ('infomax', 'extended-infomax'):
        	self.unmixing_matrix_ = infomax(Xpca, random_state=self.random_state,**self.fit_params)


        elif self.method == 'picard':
            from picard import picard
            _, W, _ = picard(Xpca.T, whiten=False,random_state=self.random_state, **self.fit_params)
            del _
            
            self.unmixing_matrix_ = W
        
        self.unmixing_matrix_ /= np.sqrt(exp_var)[None, :]  # whitening
        self.mixing_matrix_ = np.linalg.pinv(self.unmixing_matrix_)
        
        nf, ns = X.shape
        var = np.sum(self.mixing_matrix_ ** 2, axis=0) * np.sum(X**2, axis=1) / (nf*ns- 1)
        if normalize:
            var /= var.sum()

        order = var.argsort()[::-1]
        self.mixing_matrix_= self.mixing_matrix_[:, order]
        self.unmixing_matrix_ = self.unmixing_matrix_[order, :]
    def whitening(self, X):
        #pre_whitener = np.empty([len(data), 1])
        #pre_whitener[this_picks] = np.std(data[this_picks])
        whitener = np.empty([len(X), 1])
        whitener[:] = np.std(X)
        Xw = deepcopy(X)
        Xw /= whitener
        return Xw, whitener
    def transform(self, Xdata):
        """Compute sources from data (operates inplace)."""
        Xd = deepcopy(Xdata)
        if self.pca_mean_ is not None:
            Xd -= self.pca_mean_[:, None]

        # Apply first PCA
        pca_Xd = np.dot(self.pca_components_[:self.n_components], Xd)
        # Apply unmixing to low dimension PCA
        Xd = np.dot(self.unmixing_matrix_, pca_Xd)
        return Xd
    def get_sMatrix(self):
        """Get Final ICA weight matrix.
        Returns
        -------
        Matrix : array, shape (n_channels, n_components)
            The ICA weights (maps).
        """
        return np.dot(self.mixing_matrix_[:, :self.n_components].T,
                      self.pca_components_[:self.n_components]).T
    
    def get_tMatrix(self):
        return np.dot(self.unmixing_matrix_,self.pca_components_)

def infomax(data, weights=None, l_rate=None, block=None, w_change=1e-12,
            anneal_deg=60., anneal_step=0.9, extended=True, n_subgauss=1,
            kurt_size=6000, ext_blocks=1, max_iter=200, random_state=None,
            blowup=1e4, blowup_fac=0.5, n_small_angle=20, use_bias=True,
            verbose=None):
    from scipy.stats import kurtosis
    rng = check_random_state(random_state)

    # define some default parameters
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
    n_samples, n_features = data.shape
    n_features_square = n_features ** 2

    # check input parameters
    # heuristic default - may need adjustment for large or tiny data sets
    if l_rate is None:
        l_rate = 0.01 / math.log(n_features ** 2.0)

    if block is None:
        block = int(math.floor(math.sqrt(n_samples / 3.0)))

    print('Computing%sInfomax ICA' % ' Extended ' if extended else' ')

    # collect parameters
    nblock = n_samples // block
    lastt = (nblock - 1) * block + 1

    # initialize training
    if weights is None:
        weights = np.identity(n_features, dtype=np.float64)
    else:
        weights = weights.T

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
        permute = random_permutation(n_samples, rng)

        # ICA training block
        # loop across block samples
        for t in range(0, lastt, block):
            u = np.dot(data[permute[t:t + block], :], weights)
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
                        tpartact = np.dot(data[rp.astype(int), :], weights).T
                    else:
                        tpartact = np.dot(data, weights).T

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

        # here we continue after the for loop over the ICA training blocks
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
                print(
                    'step %d - lrate %5f, wchange %8.8f, angledelta %4.1f deg'
                    % (step, l_rate, change, angledelta))

            # anneal learning rate
            oldweights = weights.copy()
            if angledelta > anneal_deg:
                l_rate *= anneal_step    # anneal learning rate
                # accumulate angledelta until anneal_deg reaches l_rate
                olddelta = delta
                oldchange = change
                count_small_angle = 0  # reset count when angledelta is large
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

        # restart if weights blow up (for lowering l_rate)
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

    # prepare return values
    return weights.T

def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance.
    If seed is None, return the RandomState singleton used by np.random.
    If seed is an int, return a new RandomState instance seeded with seed.
    If seed is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (int, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)
    
def random_permutation(n_samples, random_state=None):
    """Emulate the randperm matlab function.
    It returns a vector containing a random permutation of the
    integers between 0 and n_samples-1. It returns the same random numbers
    than randperm matlab function whenever the random_state is the same
    as the matlab's random seed.
    This function is useful for comparing against matlab scripts
    which use the randperm function.
    Note: the randperm(n_samples) matlab function generates a random
    sequence between 1 and n_samples, whereas
    random_permutation(n_samples, random_state) function generates
    a random sequence between 0 and n_samples-1, that is:
    randperm(n_samples) = random_permutation(n_samples, random_state) - 1
    Parameters
    ----------
    n_samples : int
        End point of the sequence to be permuted (excluded, i.e., the end point
        is equal to n_samples-1)
    random_state : int | None
        Random seed for initializing the pseudo-random number generator.
    Returns
    -------
    randperm : ndarray, int
        Randomly permuted sequence between 0 and n-1.
    """
    rng = check_random_state(random_state)
    idx = rng.rand(n_samples)
    randperm = np.argsort(idx)
    return randperm

def check_version(library, min_version):
    r"""Check minimum library version required.
    Parameters
    ----------
    library : str
        The library name to import. Must have a ``__version__`` property.
    min_version : str
        The minimum version string. Anything that matches
        ``'(\d+ | [a-z]+ | \.)'``. Can also be empty to skip version
        check (just check for library presence).
    Returns
    -------
    ok : bool
        True if the library exists with at least the specified version.
    """
    ok = True
    try:
        library = __import__(library)
    except ImportError:
        ok = False
    else:
        if min_version:
            this_version = LooseVersion(library.__version__)
            if this_version < min_version:
                ok = False
    return ok
