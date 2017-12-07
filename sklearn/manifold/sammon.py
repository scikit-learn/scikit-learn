"""
Sammon's Mapping
"""

# author: Poom Chiarawongse <eight1911@gmail.com>
# License: BSD

# """
# He did not stop outside the door but walked quickly into the yard. His soul
# was overflowing with emotion and he felt he needed lots of room to move
# freely. Over his head was the vast vault of the sky studded with shining,
# silent stars. The still-dim Milky Way was split in two from the zenith to
# the horizon. A cool, completely still night enfolded the earth. The white
# tower and the golden domes gleamed in the sapphire sky. The gorgeous
# autumnal flowers in the flowerbeds by the building were asleep until
# morning. The silence of the earth seemed to merge with the silence of the
# sky and the mystery of the earth was one with the mystery of the stars ...
# Alyosha stood and gazed for a while; then like a blade of grass cut by a
# scythe, he fell to the ground.

# He did not know why he was hugging the earth, why he could not kiss it
# enough, why he longed to kiss it all ... He kissed it again and again,
# drenching it with his tears, vowing to love it always, always. “Water the
# earth with tears of your joy and love those tears,” a voice rang out in his
# soul. What was he weeping about? Oh, he was weeping with ecstasy, weeping,
# even, over those stars that shone down upon him from infinite distances, and
# he was “unashamed of his ecstasy.” It was as if the threads of all those
# innumerable worlds of God met in his soul and his soul was vibrating from
# its contact with “different worlds.” He craved to forgive everyone and
# everything, and to beg for forgiveness—oh, not forgiveness for just himself,
# but for everyone and everything. “Others will ask forgiveness of me too,” a
# voice rang out in his soul again. Every moment he felt clearly, almost
# physically, something real and indestructible, like the vault of the sky
# over his head, entering his soul. Something, a kind of idea, had taken over
# his soul for ever and ever. He was a weak youth when he fell on the ground
# and he rose a strong and determined fighter. He knew it. He felt it during
# that moment of rapture. And never, never thereafter would Alyosha forget
# that moment. “Someone visited my soul then,” he would say later, with firm
# faith his words.
#
# """

import numpy as np

import warnings
import random

from ..base import BaseEstimator
from ..metrics import euclidean_distances
from ..utils import check_random_state, check_array, check_symmetric

class Sammon(BaseEstimator):
    """Sammon Mapping

    Read more in the :ref:`User Guide <sammon_mapping>`.

    Parameters
    ----------
    n_components : int, optional, default: 2
        Number of dimensions in which to immerse the dissimilarities.

    l_rate : float, optional, default: 0.8
        The learning rate of Newton's gradient descent.

    decay : float, optional, default: 0.005
        The learning rate decay of Newton's gradient descent.


    base_rate : float, optional, default: 0.1
        The base learning rate. The learning rate will never drop below this.

    max_iter : int, optional, default: 300
        Maximum number of iterations of newton-raphson gradient descent.

    sensitivity : float, optional, default: 1e-5
        For converting zero entries when dividing.

    verbose : int, optional, default: 0
        Level of verbosity.

    eps : float, optional, default: 1e-3
        Relative tolerance with respect to stress at which to declare
        convergence.


    random_state : int, RandomState instance or None, optional, default: None
        The generator used to initialize the centers.  If int, random_state is
        the seed used by the random number generator; If RandomState instance,
        random_state is the random number generator; If None, the random number
        generator is the RandomState instance used by `np.random`.

    dissimilarity : 'euclidean' | 'precomputed', optional, default: 'euclidean'
        Dissimilarity measure to use:

        - 'euclidean':
            Pairwise Euclidean distances between points in the dataset.

        - 'precomputed':
            Pre-computed dissimilarities are passed directly to ``fit`` and
            ``fit_transform``.

    Attributes
    ----------
    embedding_ : array-like, shape (n_components, n_samples)
        Stores the position of the dataset in the embedding space.

    stress_ : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points).


    Reference
    ---------
    Sammon J.W. (1969) A non-linear mapping for data structure analysis.
    IEEE Transactions on Computers, 18, 401-409. pdf

    """

    def __init__(self, n_components=2, max_iter=500, l_rate=0.8, decay=0.005,
                 base_rate=0.1, verbose=0, eps=1e-3, random_state=None,
                 sensitivity=1e-5, dissimilarity="euclidean"):
        self.n_components = n_components
        self.random_state = random_state
        self.dissimilarity = dissimilarity
        self.sensitivity = sensitivity
        self.max_iter = max_iter
        self.verbose = verbose
        self.base_rate = base_rate
        self.l_rate = l_rate
        self.decay = decay
        self.eps = eps


    @property
    def _pairwise(self):
        return self.kernel == "precomputed"

    def fit(self, X, y=None, init=None):
        """
        Computes the position of the points in the embedding space

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            Input data. If ``dissimilarity=='precomputed'``, the input should
            be the dissimilarity matrix.

        y: Ignored

        init : ndarray, shape (n_samples,), optional, default: None
            Starting configuration of the embedding to initialize the SMACOF
            algorithm. By default, the algorithm is initialized with a randomly
            chosen array.
        """
        self.fit_transform(X, init=init)
        return self


    def fit_transform(self, X, y=None, init=None):
        """
        Fit the data from X, and returns the embedded coordinates

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            Input data. If ``dissimilarity=='precomputed'``, the input should
            be the dissimilarity matrix.

        y: Ignored

        init : ndarray, shape (n_samples,), optional, default: None
            Starting configuration of the embedding to initialize the SMACOF
            algorithm. By default, the algorithm is initialized with a randomly
            chosen array.
        """
        X = check_array(X)
        if X.shape[0] == X.shape[1] and self.dissimilarity != "precomputed":
            warnings.warn("The Sammon API has changed. ``fit`` "
                          "now constructs a dissimilarity matrix from data. "
                          "To use a custom dissimilarity matrix, set "
                          "``dissimilarity='precomputed'``.")

        if self.dissimilarity == "precomputed":
            self.dissimilarity_matrix_ = X
        elif self.dissimilarity == "euclidean":
            self.dissimilarity_matrix_ = euclidean_distances(X)
        else:
            raise ValueError("Proximity must be 'precomputed' or 'euclidean'."
                             " Got %s instead" % str(self.dissimilarity))

        self.embedding_, self.stress_, self.n_iter_ = sammon(
            self.dissimilarity_matrix_,
            self.n_components,
            sensitivity=self.sensitivity,
            random_state=self.random_state,
            base_rate=self.base_rate,
            max_iter=self.max_iter,
            verbose=self.verbose,
            l_rate=self.l_rate,
            decay=self.decay,
            eps=self.eps,
            init=init)

        return self.embedding_



def sammon(dissimilarity_matrix, n_components,
           init, l_rate, decay, base_rate,
           max_iter, verbose, eps, sensitivity, random_state):

    dissimilarity_matrix = check_array(dissimilarity_matrix)
    dissimilarity_matrix = check_symmetric(dissimilarity_matrix, raise_exception=True)
    random_state = check_random_state(random_state)
    n_samples = dissimilarity_matrix.shape[0]


    if init is None:
        # Randomly choose initial configuration
        X = random_state.rand(n_samples * n_components)
        X = X.reshape((n_samples, n_components))
    else:
        n_components = init.shape[1]
        if n_samples != init.shape[0]:
            raise ValueError("init matrix should be of shape (%d, %d)" %
                             (n_samples, n_components))
        X = init

    if hasattr(init, '__array__'):
        init = np.asarray(init).copy()

    pos, stress, n_iter_ = _sammon(
        dissimilarity_matrix,
        X,
        sensitivity=sensitivity,
        max_iter=max_iter,
        base_rate=base_rate,
        verbose=verbose,
        l_rate=l_rate,
        decay=decay,
        eps=eps)

    return pos, stress, n_iter_


def _sammon(dissimilarity_matrix, init, l_rate, decay, base_rate,
            max_iter, verbose, eps, sensitivity):

    n_samples = len(init)
    points = init[:]
    indices = [*range(n_samples)]
    dissimilarity_matrix = np.maximum(dissimilarity_matrix, sensitivity)

    for n_iter in range(max_iter):
        random.shuffle(indices)
        lo_dists = squareform(pdist(points))
        lo_dists = np.maximum(lo_dists, sensitivity)
        prod = dissimilarity_matrix * lo_dists
        diff = dissimilarity_matrix - lo_dists
        ratio = diff / prod
        total_delta = 0.0
        l_rate *= 1 - decay
        true_rate = l_rate + base_rate
        for a in indices:
            coord_diff = (points - points[a]).T
            grad = ratio[a] * coord_diff
            grad_prime = diff
            left = (coord_diff * coord_diff) / lo_dists[a]
            right = 1.0 + diff[a] / lo_dists[a]
            prime =  true_rate * (diff[a] * coord_diff).sum(axis=1)
            grad_prime = (diff[a] - right * left).sum(axis=1)
            delta = prime /(true_rate * np.sign(grad_prime) + grad_prime)
            points[a] += true_rate * delta
            total_delta += np.sqrt((delta**2).sum()) / n
        stress = 1 / dissimilarity_matrix.sum() * ((diff * diff) / prod).sum() / 2
        if verbose and (n_iter * verbose) % 50 == 0:
            print("iteration", n_iter, "with stress", stress, "and delta", total_delta)

        if total_delta < eps: break
    return points, stress, n_iter
