"""
Generate samples of synthetic data sets.
"""

# Authors: B. Thirion, G. Varoquaux, A. Gramfort, V. Michel, O. Grisel,
#          G. Louppe
# License: BSD 3 clause

from itertools import product
import numpy as np
from scipy import linalg

from ..utils import array2d, check_random_state


def make_classification(n_samples=100, n_features=20, n_informative=2,
                        n_redundant=2, n_repeated=0, n_classes=2,
                        n_clusters_per_class=2, weights=None, flip_y=0.01,
                        class_sep=1.0, hypercube=True, shift=0.0, scale=1.0,
                        shuffle=True, random_state=None):
    """Generate a random n-class classification problem.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    n_features : int, optional (default=20)
        The total number of features. These comprise `n_informative`
        informative features, `n_redundant` redundant features, `n_repeated`
        dupplicated features and `n_features-n_informative-n_redundant-
        n_repeated` useless features drawn at random.

    n_informative : int, optional (default=2)
        The number of informative features. Each class is composed of a number
        of gaussian clusters each located around the vertices of a hypercube
        in a subspace of dimension `n_informative`. For each cluster,
        informative features are drawn independently from  N(0, 1) and then
        randomly linearly combined in order to add covariance. The clusters
        are then placed on the vertices of the hypercube.

    n_redundant : int, optional (default=2)
        The number of redundant features. These features are generated as
        random linear combinations of the informative features.

    n_repeated : int, optional (default=2)
        The number of dupplicated features, drawn randomly from the informative
        and the redundant features.

    n_classes : int, optional (default=2)
        The number of classes (or labels) of the classification problem.

    n_clusters_per_class : int, optional (default=2)
        The number of clusters per class.

    weights : list of floats or None (default=None)
        The proportions of samples assigned to each class. If None, then
        classes are balanced. Note that if `len(weights) == n_classes - 1`,
        then the last class weight is automatically inferred.

    flip_y : float, optional (default=0.01)
        The fraction of samples whose class are randomly exchanged.

    class_sep : float, optional (default=1.0)
        The factor multiplying the hypercube dimension.

    hypercube : boolean, optional (default=True)
        If True, the clusters are put on the vertices of a hypercube. If
        False, the clusters are put on the vertices of a random polytope.

    shift : float or None, optional (default=0.0)
        Shift all features by the specified value. If None, then features
        are shifted by a random value drawn in [-class_sep, class_sep].

    scale : float or None, optional (default=1.0)
        Multiply all features by the specified value. If None, then features
        are scaled by a random value drawn in [1, 100]. Note that scaling
        happens after shifting.

    shuffle : boolean, optional (default=True)
        Shuffle the samples and the features.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, n_features]
        The generated samples.

    y : array of shape [n_samples]
        The integer labels for class membership of each sample.

    Notes
    -----
    The algorithm is adapted from Guyon [1] and was designed to generate
    the "Madelon" dataset.

    **References**:

    .. [1] I. Guyon, "Design of experiments for the NIPS 2003 variable
           selection benchmark", 2003.
    """
    generator = check_random_state(random_state)

    # Count features, clusters and samples
    assert n_informative + n_redundant + n_repeated <= n_features
    assert 2 ** n_informative >= n_classes * n_clusters_per_class
    assert weights is None or (len(weights) == n_classes or
                               len(weights) == (n_classes - 1))

    n_useless = n_features - n_informative - n_redundant - n_repeated
    n_clusters = n_classes * n_clusters_per_class

    if weights and len(weights) == (n_classes - 1):
        weights.append(1.0 - sum(weights))

    if weights is None:
        weights = [1.0 / n_classes] * n_classes
        weights[-1] = 1.0 - sum(weights[:-1])

    n_samples_per_cluster = []

    for k in xrange(n_clusters):
        n_samples_per_cluster.append(int(n_samples * weights[k % n_classes]
                                     / n_clusters_per_class))

    for i in xrange(n_samples - sum(n_samples_per_cluster)):
        n_samples_per_cluster[i % n_clusters] += 1

    # Intialize X and y
    X = np.zeros((n_samples, n_features))
    y = np.zeros(n_samples)

    # Build the polytope
    C = np.array(list(product([-class_sep, class_sep], repeat=n_informative)))

    if not hypercube:
        for k in xrange(n_clusters):
            C[k, :] *= generator.rand()

        for f in xrange(n_informative):
            C[:, f] *= generator.rand()

    generator.shuffle(C)

    # Loop over all clusters
    pos = 0
    pos_end = 0

    for k in xrange(n_clusters):
        # Number of samples in cluster k
        n_samples_k = n_samples_per_cluster[k]

        # Define the range of samples
        pos = pos_end
        pos_end = pos + n_samples_k

        # Assign labels
        y[pos:pos_end] = k % n_classes

        # Draw features at random
        X[pos:pos_end, :n_informative] = generator.randn(n_samples_k,
                                                         n_informative)

        # Multiply by a random matrix to create co-variance of the features
        A = 2 * generator.rand(n_informative, n_informative) - 1
        X[pos:pos_end, :n_informative] = np.dot(X[pos:pos_end, :n_informative],
                                                A)

        # Shift the cluster to a vertice
        X[pos:pos_end, :n_informative] += np.tile(C[k, :], (n_samples_k, 1))

    # Create redundant features
    if n_redundant > 0:
        B = 2 * generator.rand(n_informative, n_redundant) - 1
        X[:, n_informative:n_informative + n_redundant] = \
                                            np.dot(X[:, :n_informative], B)

    # Repeat some features
    if n_repeated > 0:
        n = n_informative + n_redundant
        indices = ((n - 1) * generator.rand(n_repeated) + 0.5).astype(np.int)
        X[:, n:n + n_repeated] = X[:, indices]

    # Fill useless features
    X[:, n_features - n_useless:] = generator.randn(n_samples, n_useless)

    # Randomly flip labels
    if flip_y >= 0.0:
        for i in xrange(n_samples):
            if generator.rand() < flip_y:
                y[i] = generator.randint(n_classes)

    # Randomly shift and scale
    constant_shift = shift is not None
    constant_scale = scale is not None

    for f in xrange(n_features):
        if not constant_shift:
            shift = (2 * generator.rand() - 1) * class_sep

        if not constant_scale:
            scale = 1 + 100 * generator.rand()

        X[:, f] += shift
        X[:, f] *= scale

    # Randomly permute samples and features
    if shuffle:
        indices = range(n_samples)
        generator.shuffle(indices)
        X = X[indices]
        y = y[indices]

        indices = range(n_features)
        generator.shuffle(indices)
        X[:, :] = X[:, indices]

    return X, y


def make_multilabel_classification(n_samples=100, n_features=20, n_classes=5,
                                   n_labels=2, length=50,
                                   allow_unlabeled=True, random_state=None):
    """Generate a random multilabel classification problem.

    For each sample, the generative process is:
        - pick the number of labels: n ~ Poisson(n_labels)
        - n times, choose a class c: c ~ Multinomial(theta)
        - pick the document length: k ~ Poisson(length)
        - k times, choose a word: w ~ Multinomial(theta_c)

    In the above process, rejection sampling is used to make sure that
    n is never zero or more than `n_classes`, and that the document length
    is never zero. Likewise, we reject classes which have already been chosen.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    n_features : int, optional (default=20)
        The total number of features.

    n_classes : int, optional (default=5)
        The number of classes of the classification problem.

    n_labels : int, optional (default=2)
        The average number of labels per instance. Number of labels follows
        a Poisson distribution that never takes the value 0.

    length : int, optional (default=50)
        Sum of the features (number of words if documents).

    allow_unlabeled : bool, optional (default=True)
        If ``True``, some instances might not belong to any class.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, n_features]
        The generated samples.

    Y : list of tuples
        The label sets.
    """
    generator = check_random_state(random_state)
    p_c = generator.rand(n_classes)
    p_c /= p_c.sum()
    p_w_c = generator.rand(n_features, n_classes)
    p_w_c /= np.sum(p_w_c, axis=0)

    def sample_example():
        _, n_classes = p_w_c.shape

        # pick a nonzero number of labels per document by rejection sampling
        n = n_classes + 1
        while (not allow_unlabeled and n == 0) or n > n_classes:
            n = generator.poisson(n_labels)

        # pick n classes
        y = []
        while len(y) != n:
            # pick a class with probability P(c)
            c = generator.multinomial(1, p_c).argmax()

            if not c in y:
                y.append(c)

        # pick a non-zero document length by rejection sampling
        k = 0
        while k == 0:
            k = generator.poisson(length)

        # generate a document of length k words
        x = np.zeros(n_features, dtype=int)
        for i in range(k):
            if len(y) == 0:
                # if sample does not belong to any class, generate noise word
                w = generator.randint(n_features)
            else:
                # pick a class and generate an appropriate word
                c = y[generator.randint(len(y))]
                w = generator.multinomial(1, p_w_c[:, c]).argmax()
            x[w] += 1

        return x, y

    X, Y = zip(*[sample_example() for i in range(n_samples)])
    return np.array(X, dtype=np.float64), Y


def make_regression(n_samples=100, n_features=100, n_informative=10, bias=0.0,
                    effective_rank=None, tail_strength=0.5, noise=0.0,
                    shuffle=True, coef=False, random_state=None):
    """Generate a random regression problem.

    The input set can either be well conditioned (by default) or have a low
    rank-fat tail singular profile. See the `make_low_rank_matrix` for
    more details.

    The output is generated by applying a (potentially biased) random linear
    regression model with `n_informative` nonzero regressors to the previously
    generated input and some gaussian centered noise with some adjustable
    scale.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    n_features : int, optional (default=100)
        The number of features.

    n_informative : int, optional (default=10)
        The number of informative features, i.e., the number of features used
        to build the linear model used to generate the output.

    bias : float, optional (default=0.0)
        The bias term in the underlying linear model.

    effective_rank : int or None, optional (default=None)
        if not None:
            The approximate number of singular vectors required to explain most
            of the input data by linear combinations. Using this kind of
            singular spectrum in the input allows the generator to reproduce
            the correlations often observed in practice.
        if None:
            The input set is well conditioned, centered and gaussian with
            unit variance.

    tail_strength : float between 0.0 and 1.0, optional (default=0.5)
        The relative importance of the fat noisy tail of the singular values
        profile if `effective_rank` is not None.

    noise : float, optional (default=0.0)
        The standard deviation of the gaussian noise applied to the output.

    shuffle : boolean, optional (default=True)
        Shuffle the samples and the features.

    coef : boolean, optional (default=False)
        If True, the coefficients of the underlying linear model are returned.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, n_features]
        The input samples.

    y : array of shape [n_samples]
        The output values.

    coef : array of shape [n_features], optional
        The coefficient of the underlying linear model. It is returned only if
        coef is True.
    """
    generator = check_random_state(random_state)

    if effective_rank is None:
        # Randomly generate a well conditioned input set
        X = generator.randn(n_samples, n_features)

    else:
        # Randomly generate a low rank, fat tail input set
        X = make_low_rank_matrix(n_samples=n_samples,
                                 n_features=n_features,
                                 effective_rank=effective_rank,
                                 tail_strength=tail_strength,
                                 random_state=generator)

    # Generate a ground truth model with only n_informative features being non
    # zeros (the other features are not correlated to y and should be ignored
    # by a sparsifying regularizers such as L1 or elastic net)
    ground_truth = np.zeros(n_features)
    ground_truth[:n_informative] = 100 * generator.rand(n_informative)

    y = np.dot(X, ground_truth) + bias

    # Add noise
    if noise > 0.0:
        y += generator.normal(scale=noise, size=y.shape)

    # Randomly permute samples and features
    if shuffle:
        indices = range(n_samples)
        generator.shuffle(indices)
        X = X[indices]
        y = y[indices]

        indices = range(n_features)
        generator.shuffle(indices)
        X[:, :] = X[:, indices]
        ground_truth = ground_truth[indices]

    if coef:
        return X, y, ground_truth

    else:
        return X, y


def make_blobs(n_samples=100, n_features=2, centers=3, cluster_std=1.0,
               center_box=(-10.0, 10.0), shuffle=True, random_state=None):
    """Generate isotropic Gaussian blobs for clustering.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The total number of points equally divided among clusters.

    n_features : int, optional (default=2)
        The number of features for each sample.

    centers : int or array of shape [n_centers, n_features], optional
        (default=3)
        The number of centers to generate, or the fixed center locations.

    cluster_std: float or sequence of floats, optional (default=1.0)
        The standard deviation of the clusters.

    center_box: pair of floats (min, max), optional (default=(-10.0, 10.0))
        The bounding box for each cluster center when centers are
        generated at random.

    shuffle : boolean, optional (default=True)
        Shuffle the samples.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, n_features]
        The generated samples.

    y : array of shape [n_samples]
        The integer labels for cluster membership of each sample.

    Examples
    --------
    >>> from sklearn.datasets.samples_generator import make_blobs
    >>> X, y = make_blobs(n_samples=10, centers=3, n_features=2,
    ...                   random_state=0)
    >>> X.shape
    (10, 2)
    >>> y
    array([0, 0, 1, 0, 2, 2, 2, 1, 1, 0])
    """
    generator = check_random_state(random_state)

    if isinstance(centers, int):
        centers = generator.uniform(center_box[0], center_box[1],
                                    size=(centers, n_features))
    else:
        centers = array2d(centers)
        n_features = centers.shape[1]

    X = []
    y = []

    n_centers = centers.shape[0]
    n_samples_per_center = [int(n_samples // n_centers)] * n_centers

    for i in xrange(n_samples % n_centers):
        n_samples_per_center[i] += 1

    for i, n in enumerate(n_samples_per_center):
        X.append(centers[i] + generator.normal(scale=cluster_std,
                                               size=(n, n_features)))
        y += [i] * n

    X = np.concatenate(X)
    y = np.array(y)

    if shuffle:
        indices = np.arange(n_samples)
        generator.shuffle(indices)
        X = X[indices]
        y = y[indices]

    return X, y


def make_friedman1(n_samples=100, n_features=10, noise=0.0, random_state=None):
    """Generate the "Friedman #1" regression problem

    This dataset is described in Friedman [1] and Breiman [2].

    Inputs `X` are independent features uniformly distributed on the interval
    [0, 1]. The output `y` is created according to the formula::

        y(X) = 10 * sin(pi * X[:, 0] * X[:, 1]) + 20 * (X[:, 2] - 0.5) ** 2 \
               + 10 * X[:, 3] + 5 * X[:, 4] + noise * N(0, 1).

    Out of the `n_features` features, only 5 are actually used to compute
    `y`. The remaining features are independent of `y`.

    The number of features has to be >= 5.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    n_features : int, optional (default=10)
        The number of features. Should be at least 5.

    noise : float, optional (default=0.0)
        The standard deviation of the gaussian noise applied to the output.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, n_features]
        The input samples.

    y : array of shape [n_samples]
        The output values.

    Notes
    -----
    **References**:

    .. [1] J. Friedman, "Multivariate adaptive regression splines", The Annals
           of Statistics 19 (1), pages 1-67, 1991.

    .. [2] L. Breiman, "Bagging predictors", Machine Learning 24,
           pages 123-140, 1996.
    """
    assert n_features >= 5

    generator = check_random_state(random_state)

    X = generator.rand(n_samples, n_features)
    y = 10 * np.sin(np.pi * X[:, 0] * X[:, 1]) + 20 * (X[:, 2] - 0.5) ** 2 \
        + 10 * X[:, 3] + 5 * X[:, 4] + noise * generator.randn(n_samples)

    return X, y


def make_friedman2(n_samples=100, noise=0.0, random_state=None):
    """Generate the "Friedman #2" regression problem

    This dataset is described in Friedman [1] and Breiman [2].

    Inputs `X` are 4 independent features uniformly distributed on the
    intervals::

        0 <= X[:, 0] <= 100,
        40 * pi <= X[:, 1] <= 560 * pi,
        0 <= X[:, 2] <= 1,
        1 <= X[:, 3] <= 11.

    The output `y` is created according to the formula::

        y(X) = (X[:, 0] ** 2 \
                   + (X[:, 1] * X[:, 2] \
                         - 1 / (X[:, 1] * X[:, 3])) ** 2) ** 0.5 \
               + noise * N(0, 1).

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    noise : float, optional (default=0.0)
        The standard deviation of the gaussian noise applied to the output.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, 4]
        The input samples.

    y : array of shape [n_samples]
        The output values.

    Notes
    -----
    **References**:

    .. [1] J. Friedman, "Multivariate adaptive regression splines", The Annals
           of Statistics 19 (1), pages 1-67, 1991.

    .. [2] L. Breiman, "Bagging predictors", Machine Learning 24,
           pages 123-140, 1996.
    """
    generator = check_random_state(random_state)

    X = generator.rand(n_samples, 4)
    X[:, 0] *= 100
    X[:, 1] *= 520 * np.pi
    X[:, 1] += 40 * np.pi
    X[:, 3] *= 10
    X[:, 3] += 1

    y = (X[:, 0] ** 2
            + (X[:, 1] * X[:, 2] - 1 / (X[:, 1] * X[:, 3])) ** 2) ** 0.5 \
        + noise * generator.randn(n_samples)

    return X, y


def make_friedman3(n_samples=100, noise=0.0, random_state=None):
    """Generate the "Friedman #3" regression problem

    This dataset is described in Friedman [1] and Breiman [2].

    Inputs `X` are 4 independent features uniformly distributed on the
    intervals::

        0 <= X[:, 0] <= 100,
        40 * pi <= X[:, 1] <= 560 * pi,
        0 <= X[:, 2] <= 1,
        1 <= X[:, 3] <= 11.

    The output `y` is created according to the formula::

        y(X) = arctan((X[:, 1] * X[:, 2] \
                          - 1 / (X[:, 1] * X[:, 3])) \
                      / X[:, 0]) \
               + noise * N(0, 1).

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    noise : float, optional (default=0.0)
        The standard deviation of the gaussian noise applied to the output.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, 4]
        The input samples.

    y : array of shape [n_samples]
        The output values.

    Notes
    -----
    **References**:

    .. [1] J. Friedman, "Multivariate adaptive regression splines", The Annals
           of Statistics 19 (1), pages 1-67, 1991.

    .. [2] L. Breiman, "Bagging predictors", Machine Learning 24,
           pages 123-140, 1996.
    """
    generator = check_random_state(random_state)

    X = generator.rand(n_samples, 4)
    X[:, 0] *= 100
    X[:, 1] *= 520 * np.pi
    X[:, 1] += 40 * np.pi
    X[:, 3] *= 10
    X[:, 3] += 1

    y = np.arctan((X[:, 1] * X[:, 2] - 1 / (X[:, 1] * X[:, 3])) / X[:, 0]) \
        + noise * generator.randn(n_samples)

    return X, y


def make_low_rank_matrix(n_samples=100, n_features=100, effective_rank=10,
                         tail_strength=0.5, random_state=None):
    """Generate a mostly low rank matrix with bell-shaped singular values

    Most of the variance can be explained by a bell-shaped curve of width
    effective_rank: the low rank part of the singular values profile is::

        (1 - tail_strength) * exp(-1.0 * (i / effective_rank) ** 2)

    The remaining singular values' tail is fat, decreasing as::

        tail_strength * exp(-0.1 * i / effective_rank).

    The low rank part of the profile can be considered the structured
    signal part of the data while the tail can be considered the noisy
    part of the data that cannot be summarized by a low number of linear
    components (singular vectors).

    This kind of singular profiles is often seen in practice, for instance:
     - gray level pictures of faces
     - TF-IDF vectors of text documents crawled from the web

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    n_features : int, optional (default=100)
        The number of features.

    effective_rank : int, optional (default=10)
        The approximate number of singular vectors required to explain most of
        the data by linear combinations.

    tail_strength : float between 0.0 and 1.0, optional (default=0.5)
        The relative importance of the fat noisy tail of the singular values
        profile.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, n_features]
        The matrix.
    """
    generator = check_random_state(random_state)
    n = min(n_samples, n_features)

    # Random (ortho normal) vectors
    from ..utils.fixes import qr_economic
    u, _ = qr_economic(generator.randn(n_samples, n))
    v, _ = qr_economic(generator.randn(n_features, n))

    # Index of the singular values
    singular_ind = np.arange(n, dtype=np.float64)

    # Build the singular profile by assembling signal and noise components
    low_rank = (1 - tail_strength) * \
               np.exp(-1.0 * (singular_ind / effective_rank) ** 2)
    tail = tail_strength * np.exp(-0.1 * singular_ind / effective_rank)
    s = np.identity(n) * (low_rank + tail)

    return np.dot(np.dot(u, s), v.T)


def make_sparse_coded_signal(n_samples, n_components, n_features,
                             n_nonzero_coefs, random_state=None):
    """Generate a signal as a sparse combination of dictionary elements.

    Returns a matrix Y = DX, such as D is (n_features, n_components),
    X is (n_components, n_samples) and each column of X has exactly
    n_nonzero_coefs non-zero elements.

    Parameters
    ----------
    n_samples : int
        number of samples to generate

    n_components:  int,
        number of components in the dictionary

    n_features : int
        number of features of the dataset to generate

    n_nonzero_coefs : int
        number of active (non-zero) coefficients in each sample

    random_state: int or RandomState instance, optional (default=None)
        seed used by the pseudo random number generator

    Returns
    -------
    data: array of shape [n_features, n_samples]
        The encoded signal (Y).

    dictionary: array of shape [n_features, n_components]
        The dictionary with normalized components (D).

    code: array of shape [n_components, n_samples]
        The sparse code such that each column of this matrix has exactly
        n_nonzero_coefs non-zero items (X).

    """
    generator = check_random_state(random_state)

    # generate dictionary
    D = generator.randn(n_features, n_components)
    D /= np.sqrt(np.sum((D ** 2), axis=0))

    # generate code
    X = np.zeros((n_components, n_samples))
    for i in xrange(n_samples):
        idx = np.arange(n_components)
        generator.shuffle(idx)
        idx = idx[:n_nonzero_coefs]
        X[idx, i] = generator.randn(n_nonzero_coefs)

    # encode signal
    Y = np.dot(D, X)

    return map(np.squeeze, (Y, D, X))


def make_sparse_uncorrelated(n_samples=100, n_features=10, random_state=None):
    """Generate a random regression problem with sparse uncorrelated design

    This dataset is described in Celeux et al [1]. as::

        X ~ N(0, 1)
        y(X) = X[:, 0] + 2 * X[:, 1] - 2 * X[:, 2] - 1.5 * X[:, 3]

    Only the first 4 features are informative. The remaining features are
    useless.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    n_features : int, optional (default=10)
        The number of features.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, n_features]
        The input samples.

    y : array of shape [n_samples]
        The output values.

    Notes
    -----
    **References**:

    .. [1] G. Celeux, M. El Anbari, J.-M. Marin, C. P. Robert,
           "Regularization in regression: comparing Bayesian and frequentist
           methods in a poorly informative situation", 2009.
    """
    generator = check_random_state(random_state)

    X = generator.normal(loc=0, scale=1, size=(n_samples, n_features))
    y = generator.normal(loc=(X[:, 0] +
                              2 * X[:, 1] -
                              2 * X[:, 2] -
                              1.5 * X[:, 3]), scale=np.ones(n_samples))

    return X, y


def make_spd_matrix(n_dim, random_state=None):
    """Generate a random symmetric, positive-definite matrix.

    Parameters
    ----------
    n_dim : int
        The matrix dimension.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_dim, n_dim]
        The random symmetric, positive-definite matrix.
    """
    generator = check_random_state(random_state)

    A = generator.rand(n_dim, n_dim)
    U, s, V = linalg.svd(np.dot(A.T, A))
    X = np.dot(np.dot(U, 1.0 + np.diag(generator.rand(n_dim))), V)

    return X


def make_sparse_spd_matrix(dim=1, alpha=0.95, norm_diag=False,
                           smallest_coef=.1, largest_coef=.9,
                           random_state=None):
    """Generate a sparse symetric definite positive matrix.

    Parameters
    ----------
    dim: integer, optional (default=1)
        The size of the random  (matrix to generate.

    alpha: float between 0 and 1, optional (default=0.95)
        The probability that a coefficient is non zero (see notes).

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    prec: array of shape = [dim, dim]

    Notes
    -----
    The sparsity is actually imposed on the cholesky factor of the matrix.
    Thus alpha does not translate directly into the filling fraction of
    the matrix itself.
    """
    random_state = check_random_state(random_state)

    chol = -np.eye(dim)
    aux = random_state.rand(dim, dim)
    aux[aux < alpha] = 0
    aux[aux > alpha] = (smallest_coef
                        + (largest_coef - smallest_coef)
                          * random_state.rand(np.sum(aux > alpha)))
    aux = np.tril(aux, k=-1)

    # Permute the lines: we don't want to have assymetries in the final
    # SPD matrix
    permutation = random_state.permutation(dim)
    aux = aux[permutation].T[permutation]
    chol += aux
    prec = np.dot(chol.T, chol)

    if norm_diag:
        d = np.diag(prec)
        d = 1. / np.sqrt(d)
        prec *= d
        prec *= d[:, np.newaxis]

    return prec


def make_swiss_roll(n_samples=100, noise=0.0, random_state=None):
    """Generate a swiss roll dataset.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of sample points on the S curve.

    noise : float, optional (default=0.0)
        The standard deviation of the gaussian noise.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, 3]
        The points.

    t : array of shape [n_samples]
        The univariate position of the sample according to the main dimension
        of the points in the manifold.

    Notes
    -----
    The algorithm is from Marsland [1].

    **References**:

    .. [1] S. Marsland, "Machine Learning: An Algorithmic Perpsective",
           Chapter 10, 2009.
           http://www-ist.massey.ac.nz/smarsland/Code/10/lle.py
    """
    generator = check_random_state(random_state)

    t = 1.5 * np.pi * (1 + 2 * generator.rand(1, n_samples))
    x = t * np.cos(t)
    y = 21 * generator.rand(1, n_samples)
    z = t * np.sin(t)

    X = np.concatenate((x, y, z))
    X += noise * generator.randn(3, n_samples)
    X = X.T
    t = np.squeeze(t)

    return X, t


def make_s_curve(n_samples=100, noise=0.0, random_state=None):
    """Generate an S curve dataset.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of sample points on the S curve.

    noise : float, optional (default=0.0)
        The standard deviation of the gaussian noise.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, 3]
        The points.

    t : array of shape [n_samples]
        The univariate position of the sample according to the main dimension
        of the points in the manifold.
    """
    generator = check_random_state(random_state)

    t = 3 * np.pi * (generator.rand(1, n_samples) - 0.5)
    x = np.sin(t)
    y = 2.0 * generator.rand(1, n_samples)
    z = np.sign(t) * (np.cos(t) - 1)

    X = np.concatenate((x, y, z))
    X += noise * generator.randn(3, n_samples)
    X = X.T
    t = np.squeeze(t)

    return X, t
