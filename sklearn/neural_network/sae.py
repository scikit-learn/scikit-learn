import numpy as np
from scipy.optimize import fmin_l_bfgs_b
from scipy.linalg import norm
from ..utils import array2d,check_random_state
from sklearn.base import BaseEstimator, TransformerMixin

def sigmoid(x):
    """
    Implements the sigmoid function.

    Parameters
    ----------
    x: array-like, shape (M, N)

    Returns
    -------
    x_new: array-like, shape (M, N)
    """
    return 1. / (1. + np.exp(np.clip(-x,-30,30)))
    

class SAE(BaseEstimator, TransformerMixin):
    """
    Sparse Autoencoder (SAE)

    A Sparse Autoencoder with one hidden layer. Parameters are trained using 
    Limited-Memory  Broyden-Fletcher-Goldfarb-Shanno (L-BFGS)
    Parameters
    ----------
    n_hidden : int
        Number of hidden neurons
    lr : float, optional
        Learning rate to use during learning. It is *highly* recommended
        to tune this hyper-parameter. Possible values are 10**[0., -3.].
    beta : float, optional
        Weight of sparsity penalty term   
    sparsityParam : float, optional
        Desired average activation of the hidden units
    batch_size : int, optional
        Number of examples per minibatch.
    n_iter : int, optional
        Number of iterations/sweeps over the training dataset to perform
        during training.
    verbose: bool, optional
        When True (False by default) the method outputs the progress
        of learning after each iteration.
    random_state : integer or numpy.RandomState, optional
        A random number generator instance to define the state of the
        random permutations generator. If an integer is given, it fixes the
        seed. Defaults to the global numpy random number generator.

    Attributes
    ----------
    self.coef_hidden_ : array-like, shape (n_hidden, n_visible)
        Weight matrix, where n_visible in the number of visible
        units and n_hidden is the number of hidden units.
    self.coef_output_  : array-like, shape (n_visible, n_hidden)      
        Weight matrix, where n_visible in the number of visible
        units and n_hidden is the number of hidden units.    
    intercept_hidden_ : array-like, shape (n_hidden,), optional
        Biases of the hidden units
    intercept_visible_ : array-like, shape (n_visible,), optional
        Biases of the visible units

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.neural_network import SAE
    >>> X = np.array([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1]])
    >>> model = SAE(n_hidden=10)
    >>> model.fit(X)
    SAE(batch_size=10000, beta=3, lr=0.0001, n_hidden=10, n_iter=20,
  random_state=None, sparsityParam=0.01, verbose=False)

    References
    ----------

    [1] Ngiam, Jiquan, et al. "On optimization methods for deep learning."
        Proceedings of the 28th International Conference on Machine Learning (ICML-11). 2011.
        http://ai.stanford.edu/~quocle/LeNgiCoaLahProNg11.pdf
    """
    def __init__(self, n_hidden= 25, lr = 0.0001, beta = 3, sparsityParam = 0.01, batch_size = 10000,n_iter=20, verbose=False,random_state=None):
        self.n_hidden = n_hidden
        self.lr = lr
        self.beta = beta
        self.sparsityParam = sparsityParam
        self.batch_size = batch_size
        self.n_iter=n_iter
        self.verbose = verbose
        self.random_state=random_state
    def transform(self,X):
        """
        Computes the extracted features.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)

        Returns
        -------
        h: array-like, shape (n_samples, n_components)
        """
        return sigmoid(np.dot(X,self.coef_hidden_.T) + self.intercept_hidden_)
    def fit_transform(self,X):
        """
        Fit the model to the data X and transform it.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.
        """
        self.fit(X)
        return self.transform(X)
    def fit(self,X):
        """
        Fit the model to the data X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self
        """
        X = array2d(X)
        dtype = np.float32 if X.dtype.itemsize == 4 else np.float64
        rng = check_random_state(self.random_state)
        initial_theta=self._initParams(X.shape[1])
        inds = np.arange(X.shape[0])
        rng.shuffle(inds)
        n_batches = int(np.ceil(len(inds) / float(self.batch_size)))
        #Having input X transposed improves the performance substantially
        X=X.T
        for minibatch in xrange(n_batches):
                self._fit(X[:,inds[minibatch::n_batches]],initial_theta)
        return self
    def _fit(self,X,initial_theta):
        """
        Train the parameters using Limited-Memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) [1].

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.
        initial_theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2), 1)
            Contains concatenated flattened weights  that represent the parameters "W1, W2, b1, b2"

        """
        n_samples=X.shape[1]
        n_visible=X.shape[0]
        options = { 'maxfun':  self.n_iter,'disp':self.verbose}
        #fmin_l_bfgs_b is a very powerful and efficient optimizer
        optTheta, cost, d = \
            fmin_l_bfgs_b(lambda t: self._cost(t,X, n_visible,n_samples),
                    initial_theta, **options)
        #store the optimal weights and biases in coefficient and intercept attributes
        self._unpack(optTheta,n_visible)
         
    def _initParams(self,n_visible):
        """
        Initialize weight and bias parameters

        Parameters
        ----------
        n_visible: int
            Number of features (visible nodes).
        
        Returns
        -------
        theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2), 1)
        """
        #TODO: use 'rng' to set randomness seed
        # choose uniformly in the interval [-r, r]
        r  = np.sqrt(6) / np.sqrt(self.n_hidden+n_visible+1);   
        W1 = np.random.rand(self.n_hidden, n_visible) * 2 * r - r;
        W2 = np.random.rand(n_visible, self.n_hidden) * 2 * r - r;
        b1 = np.zeros(self.n_hidden);
        b2 = np.zeros(n_visible);
        return np.hstack((W1.ravel(), W2.ravel(), b1.ravel(), b2.ravel()))
    def _unpack(self,theta,n_visible):
          """
          Extract the coefficients and intercepts (W1,W2,b1,b2) from theta

          Parameters
          ----------
          theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2), 1)
            Contains concatenated flattened weights  that represent the parameters "W1, W2, b1, b2"
          n_visible: int
            Number of features (visible nodes).
          """
          N = self.n_hidden * n_visible
          self.coef_hidden_ = np.reshape(theta[:N],
              (self.n_hidden, n_visible))
          self.coef_output_ = np.reshape(theta[N:2*N],
              (n_visible, self.n_hidden))
          self.intercept_hidden_ = theta[2*N:2*N + self.n_hidden]
          self.intercept_output_ = theta[2*N + self.n_hidden:]
    def _cost(self,theta,X,n_visible,n_samples):
          """
          Computes the sparse autoencoder cost  function ``Jsparse(W,b)``
          and the corresponding derivatives of Jsparse with respect to the 
          different parameters given in the initialization [1]

          Parameters
          ----------
          theta: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2))
            Contains concatenated flattened weights  that represent the parameters "W1, W2, b1, b2"
          X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.
          n_visible: int
            Number of features (visible nodes).
          n_samples: int
            Number of samples
            
         Returns
         -------
         cost: float
         grad: array-like, shape (size(W1)*size(W2)*size(b1)*size(b2))
         
         References
         -------
         [1] http://ufldl.stanford.edu/wiki/index.php/Autoencoders_and_Sparsity
          """
          # Unpack theta to extract W1, W2, b1, b2
          self._unpack(theta,n_visible)
         # Forward propagate
          a2 = sigmoid(np.dot(self.coef_hidden_, X) + self.intercept_hidden_[:,np.newaxis])
          a3 = sigmoid(np.dot(self.coef_output_, a2) + self.intercept_output_[:,np.newaxis])
          # Get average activation of hidden neurons
          sparsityParamHat = np.mean(a2, 1)
          sparsity_delta  = self.beta * ((1 - self.sparsityParam) / (1 - sparsityParamHat) - self.sparsityParam / sparsityParamHat)
          # Backward propagate
          diff = X - a3
          delta3 = -diff * a3 * (1 - a3)
          delta2 = ((np.dot(self.coef_output_.T, delta3) +  sparsity_delta[:,np.newaxis]))* a2 * (1 - a2)
          # Get cost and gradient
          cost = np.trace(np.dot(diff, diff.T)) / (2 * n_samples)
          W1grad = np.dot(delta2, X.T) / n_samples
          W2grad = np.dot(delta3 , a2.T) / n_samples
          b1grad = np.sum(delta2, 1) / n_samples
          b2grad = np.sum(delta3, 1) / n_samples
          # Add regularization term to cost and gradient
          cost = cost + (self.lr / 2) * (norm(self.coef_hidden_, 'fro')**2 + norm(self.coef_output_, 'fro')**2)
          W1grad = W1grad + self.lr * self.coef_hidden_
          W2grad = W2grad + self.lr * self.coef_output_
          # Add sparsity term to the cost
          sparseCost = np.sum(self.sparsityParam * np.log(self.sparsityParam / sparsityParamHat) + (1 - self.sparsityParam) *
              np.log((1-self.sparsityParam)/(1-sparsityParamHat)))
          cost = cost + self.beta * sparseCost
          # Convert grad to vector form (This is necessary for the fmin optimizer)
          grad = np.hstack((W1grad.ravel(), W2grad.ravel(),
                  b1grad.ravel(), b2grad.ravel()))
          return cost, grad



