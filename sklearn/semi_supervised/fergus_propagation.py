import numpy as np
import scipy.linalg as LA
import scipy.sparse
import sklearn.utils.arpack as SLA
from sklearn.base import ClassifierMixin
from sklearn.base import BaseEstimator
import sklearn.metrics.pairwise as pairwise
import sklearn.mixture as mixture

class FergusPropagation(BaseEstimator, ClassifierMixin):
    '''
    The Fergus et al 2009 eigenfunction propagation classifier.

    Parameters
    ----------
    kernel : {'rbf', 'img', 'knn'}
        String identifier for kernel function to use.
        rbf: Creates a fully-connected dense RBF kernel.
        img: Creates a sparse RBF kernel where each pixel is connected only
            to its 8-neighborhood; all others are 0. ONLY use for image data!
            NOTE: Currently only support .png images.
        knn: Creates RBF kernel with connections to closest n_neighbors points,
            measured in euclidean distance.
    k : integer > 0
        Number of eigenvectors to use (default: # of labels, not including 
        unlabeled data) in eigen-decomposition.
    gamma : float > 0
        Parameter for RBF kernel (modes: rbf, img) (default: 20).
    n_neighbors : int > 0
        Parameter for the KNN kernel (modes: knn).
    lagrangian : float > 0
        Lagrangian multiplier to constrain the regularization penalty (default: 10).
    img_dims : array-like, shape = [width, height]
        Tuple specifying the width and height (or COLS and ROWS, respectively)
        of the input image. Used only with a 'img' kernel type.

    Examples
    --------
    TODO

    References
    ----------
    Rob Fergus, Yair Weiss, Antonio Torralba. Semi-Supervised Learning in 
    Gigantic Image Collections (2009).
    http://eprints.pascal-network.org/archive/00005636/01/ssl-1.pdf
    '''

    def __init__(self, kernel = 'rbf', k = -1, gamma = 20, n_neighbors = 7, lagrangian = 10, img_dims = (-1, -1)):
        # This doesn't iterate at all, so the parameters are very different
        # from the BaseLabelPropagation parameters.
        self.kernel = kernel
        self.k = k
        self.gamma = gamma
        self.n_neighbors = n_neighbors
        self.lagrangian = lagrangian
        self.img_dims = img_dims

    def fit(self, X, y):
        '''
        Fit a semi-supervised eigenfunction label propagation model.

        All input data is provided in matrix X (labeled and unlabeled), and
        corresponding label vector y with dedicated marker value (-1) for 
        unlabeled samples.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Feature vectors for data.
        y : array-like, shape = [n_samples]
            Label array, unlabeled points are marked as -1. All unlabeled
            samples will be assigned labels.

        Returns
        -------
        self : An instance of self.
        '''
        self.X_ = X

        # Construct the graph laplacian from the graph matrix.
        W = self._get_kernel(self.X_)
        L = self._normalized_graph_laplacian(W)

        # Perform the eigen-decomposition.
        vals = None
        vects = None
        classes = np.sort(np.unique(y))
        if classes[0] == -1:
            classes = np.delete(classes, 0) # remove the -1 from this list
            if self.k == -1:
                self.k = np.size(classes)

        if scipy.sparse.isspmatrix(L):
            vals, vects = SLA.eigsh(L, k = self.k, which = 'LM', maxiter = 1000)
        else:
            M = L.shape[0]
            vals, vects = LA.eigh(L, eigvals = (M - self.k, M - 1))
        self.u_ = vals

        # Construct some matrices.
        # U: k eigenvectors corresponding to smallest eigenvalues. (n_samples by k)
        # S: Diagonal matrix of k smallest eigenvalues. (k by k)
        # V: Diagonal matrix of LaGrange multipliers for labeled data, 0 for unlabeled. (n_samples by n_samples)
        self.U_ = np.real(vects)
        S = np.diag(self.u_)
        V = np.diag(np.zeros(np.size(y)))
        labeled = np.where(y != -1)
        V[labeled, labeled] = self.lagrangian

        # Solve for alpha and use it to compute the eigenfunctions, f.
        A = S + np.dot(np.dot(self.U_.T, V), self.U_)
        b = np.dot(np.dot(self.U_.T, V), y)
        alpha = LA.solve(A, b)
        f = np.dot(self.U_, alpha)

        # Set up a GMM to assign the hard labels from the eigenfunctions.
        g = mixture.GMM(n_components = np.size(classes))
        g.fit(f)
        self.labels_ = g.predict(f)
        means = np.argsort(g.means_.flatten())
        for i in range(0, np.size(self.labels_)):
            self.labels_[i] = np.where(means == self.labels_[i])[0][0]

        # Done!
        return self

	def predict(self, X):
		'''
        Performs inductive inference across the model.

        Parameters
        ----------
        X : array_like, shape = [n_samples, n_features]

        Returns
        -------
        y : array_like, shape = [n_samples]
            Predictions for input data
        '''
        probas = self.predict_proba(X)
        return self.labels_[np.argmax(probas, axis=1)].ravel()

	def predict_proba(self, X):
		'''
        Predict probability for each possible outcome.

        Compute the probability estimates for each single sample in X
        and each possible outcome seen during training (categorical
        distribution).

        Parameters
        ----------
        X : array_like, shape = [n_samples, n_features]

        Returns
        -------
        probabilities : array, shape = [n_samples, n_classes]
            Normalized probability distributions across
            class labels
        '''
        X_2d = None
        if scipy.sparse.isspmatrix(X):
            X_2d = X
        else:
            X_2d = np.atleast_2d(X)
        weight_matrices = self._get_kernel(self.X_, X_2d)
        if self.kernel == 'knn':
            probabilities = []
            for weight_matrix in weight_matrices:
                ine = np.sum(self.label_distributions_[weight_matrix], axis=0)
                probabilities.append(ine)
            probabilities = np.array(probabilities)
        else:
            weight_matrices = weight_matrices.T
            probabilities = np.dot(weight_matrices, self.label_distributions_)
        normalizer = np.atleast_2d(np.sum(probabilities, axis=1)).T
        probabilities /= normalizer
        return probabilities

    def _get_kernel(self, X, y = None):
        if self.kernel == "rbf":
            if y is None:
                return pairwise.rbf_kernel(X, X, gamma = self.gamma)
            else:
                return pairwise.rbf_kernel(X, y, gamma = self.gamma)
        elif self.kernel == "img":
            return self._img_rbf_kernel(X)
        elif self.kernel == "knn":
            if self.nn_fit is None:
                self.nn_fit = NearestNeighbors(self.n_neighbors).fit(X)
            if y is None:
                return self.nn_fit.kneighbors_graph(self.nn_fit._fit_X,
                        self.n_neighbors, mode = 'connectivity')
            else:
                return self.nn_fit.kneighbors(y, return_distance = False)
        else:
            raise ValueError("%s is not a valid kernel. Only rbf and knn"
                             " are supported at this time" % self.kernel)

    def _img_rbf_kernel(self, X):
        A = np.diag(np.ones(X.shape[0]))
        index = 0
        for i in range(0, self.img_dims[1]):
            for j in range(0, self.img_dims[0]):

                if i - 1 >= 0:
                    # Directly above.
                    other = index - self.img_dims[0]
                    rbf = pairwise.rbf_kernel(X[index], X[other], gamma = self.gamma)
                    A[index, other] = rbf
                    A[other, index] = rbf

                if i + 1 < self.img_dims[1]:
                    # Directly below.
                    other = index + self.img_dims[0]
                    rbf = pairwise.rbf_kernel(X[index], X[other], gamma = self.gamma)
                    A[index, other] = rbf
                    A[other, index] = rbf

                if j - 1 >= 0:
                    # Directly to the left.
                    other = index - 1
                    rbf = pairwise.rbf_kernel(X[index], X[other], gamma = self.gamma)
                    A[index, other] = rbf
                    A[other, index] = rbf

                if j + 1 < self.img_dims[0]:
                    # Directly to the right.
                    other = index + 1
                    rbf = pairwise.rbf_kernel(X[index], X[other], gamma = self.gamma)
                    A[index, other] = rbf
                    A[other, index] = rbf

                if i - 1 >= 0 and j - 1 >= 0:
                    # Upper left corner.
                    other = index - self.img_dims[0] - 1
                    rbf = pairwise.rbf_kernel(X[index], X[other], gamma = self.gamma)
                    A[index, other] = rbf
                    A[other, index] = rbf

                if i - 1 >= 0 and j + 1 < self.img_dims[0]:
                    # Upper right corner.
                    other = index - self.img_dims[0] + 1
                    rbf = pairwise.rbf_kernel(X[index], X[other], gamma = self.gamma)
                    A[index, other] = rbf
                    A[other, index] = rbf

                if i + 1 < self.img_dims[1] and j - 1 >= 0:
                    # Lower left corner.
                    other = index + self.img_dims[0] - 1
                    rbf = pairwise.rbf_kernel(X[index], X[other], gamma = self.gamma)
                    A[index, other] = rbf
                    A[other, index] = rbf

                if i + 1 < self.img_dims[1] and j + 1 < self.img_dims[0]:
                    # Lower right corner.
                    other = index + self.img_dims[0] + 1
                    rbf = pairwise.rbf_kernel(X[index], X[other], gamma = self.gamma)
                    A[index, other] = rbf
                    A[other, index] = rbf

                index += 1
        return scipy.sparse.csc_matrix(A)

    def _normalized_graph_laplacian(self, A):
        '''
        Calculates the normalized graph laplacian, as the sklearn utility 
        graph_laplacian(normed = True) does not seem to do this.

        Parameters
        ----------
        A : array, shape (n, n)
            Symmetric affinity matrix.

        Returns
        -------
        L : array, shape (n, n)
            Normalized graph laplacian.
        '''
        L = None
        if scipy.sparse.isspmatrix(A):

            # SKLEARN
            n_nodes = A.shape[0]
            if not A.format == 'coo':
                L = A.tocoo()
            else:
                L = A.copy()
            w = np.sqrt(np.asarray(L.sum(axis = 0)).squeeze())
            L.data /= w[L.row]
            L.data /= w[L.col]
        else:
            L = A.copy()
            d = np.sqrt(L.sum(axis = 0))
            L /= d
            L /= d[:, np.newaxis]

            # The last line effectively transposes the vector d, swapping
            # the rows and columns of the division operation. Pretty clever,
            # scikit-learn!
        return L
