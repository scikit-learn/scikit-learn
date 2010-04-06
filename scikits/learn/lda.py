# -*- coding: utf-8 -*-
import numpy as np
import exceptions
import scipy.linalg as linalg

class LDA(object):
    """
    Linear Discriminant Analysis (LDA)

    Parameters
    ----------
    X : array-like, shape = [nsamples, nfeatures]
        Training vector, where nsamples in the number of samples and
        nfeatures is the number of features.
    Y : array, shape = [nsamples]
        Target vector relative to X

    priors : array, optional, shape = [n_classes]
        ### TODO

    use_svd : bool, optional
         Specify if the SVD from scipy should be used.

    Attributes
    ----------
    `means_` : array-like, shape = [n_classes]
        Support vectors

    Methods
    -------
    fit(X, Y) : self
        Fit the model

    predict(X) : array
        Predict using the model.

    Examples
    --------
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> Y = np.array([1, 1, 1, 2, 2, 2])
    >>> clf = LDA()
    >>> clf.fit(X, Y)    #doctest: +ELLIPSIS
    <scikits.learn.lda.LDA object at 0x...>
    >>> print clf.predict([[-0.8, -1]])
    [1]

    See also
    --------
    QDA

    """
    def __init__(self, priors = None, use_svd = True):
        #use_svd : if True, use linalg.svd alse use computational trick with covariance matrix
        if not priors is None:
            self.priors = np.asarray(priors)
        else: self.priors = None
        self.use_svd = use_svd

    def fit(self, X, y, tol = 1.0e-4):
        if X.ndim!=2:
            raise exceptions.ValueError('X must be 2D array')
        n_samples = X.shape[0]
        n_features = X.shape[1]
        classes = np.unique(y)
        n_classes = classes.shape[0]
        if n_classes < 2:
            raise exceptions.ValueError('Y has Less than 2 classes')
        classes_indexes = [(y == lev).ravel() for lev in classes]
        if self.priors is None:
            counts = np.array([float(np.sum(group_indexes)) \
                for group_indexes in classes_indexes])
            self.priors = counts / n_samples
        # Group means n_classes*nfeats matrix
        means = []
        Xc = []
        for group_indexes in classes_indexes:
            Xg = X[group_indexes,:]
            meang = Xg.mean(0)
            means.append(meang)
            # centered group data
            Xgc = Xg - meang
            Xc.append(Xgc)
        means = np.asarray(means)
        Xc = np.concatenate(Xc,0)
        # ----------------------------
        # 1) within (univariate) scaling by with classes std-dev
        scaling = np.diag(1 / Xc.std(0))
        fac = float(1) / (n_samples - n_classes)
        # ----------------------------
        # 2) Within variance scaling
        X = np.sqrt(fac)*np.dot(Xc, scaling)
        # SVD of centered (within)scaled data
        if self.use_svd == True:
            U, S, V = linalg.svd(X, full_matrices=0)
        else:
            S, V = self.svd(X)

        rank = np.sum(S > tol)
        if rank < n_features: print "Warning variables are collinear"
        # Scaling of within covariance is: V' 1/S
        scaling = np.dot(np.dot(scaling, V.T[:,:rank]), np.diag(1 / S[:rank]))
        ## ----------------------------
        ## 3) Between variance scaling
        # Overall mean
        xbar = np.dot(self.priors, means)
        # Scale weighted centers
        X = np.dot(np.dot(np.diag(np.sqrt((n_samples * self.priors)*fac)),
                          (means - xbar)),
                   scaling)
        # Centers are living in a space with n_classes-1 dim (maximum)
        # Use svd to find projection in the space spamed by the (n_classes) centers
        if self.use_svd == True:
            U, S, V = linalg.svd(X,full_matrices=0)
        else:
            S, V = self.svd(X)

        rank = np.sum(S > tol*S[0])
        # compose the scalings
        scaling = np.dot(scaling,V.T[:,:rank])
        self.scaling = scaling
        self.means_ = means
        self.xbar = xbar
        self.classes = classes
        return self

    def svd(self, X):
        #computational trick to compute svd. U, S, V=linalg.svd(X)
        K = np.dot(X.T, X)
        S, V = linalg.eigh(K)
        S = np.sqrt(np.maximum(S, 1e-30))
        S_sort = -np.sort(-S)[:X.shape[0]]
        S_argsort = np.argsort(-S).tolist()
        V = V.T[S_argsort, :]
        V = V[:X.shape[0], :]
        return S_sort, V

    def predict(self, X, posterior=False):
        #Ensure X is an array
        X = np.asarray(X)
        scaling = self.scaling
        # Remove overall mean (center) and scale
        # a) data
        X = np.dot(X - self.xbar, scaling)
        # b) centers
        dm = np.dot(self.means_ - self.xbar, scaling)
        # for each class k, compute the linear discrinant function(p. 87 Hastie)
        # of sphered (scaled data)
        dist = 0.5*np.sum(dm**2, 1) - np.log(self.priors) - np.dot(X,dm.T)
        self.dist = dist
        # take exp of min dist
        dist = np.exp(-dist + dist.min(1).reshape(X.shape[0],1))
        # normalize by p(x)=sum_k p(x|k)
        self.posteriors = dist / dist.sum(1).reshape(X.shape[0],1)
        # classify according to the maximun a posteriori
        y_pred = self.classes[self.posteriors.argmax(1)]
        if posterior is True:
            return y_pred, self.posteriors
        return y_pred
