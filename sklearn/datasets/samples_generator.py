"""
Generate samples of synthetic data sets.
"""

# Authors: B. Thirion, G. Varoquaux, A. Gramfort, V. Michel, O. Grisel,
#          G. Louppe, J. Nothman
# License: BSD 3 clause

import numbers
import array
import numpy as np
from scipy import linalg
import scipy.sparse as sp

from ..preprocessing import MultiLabelBinarizer
from ..utils import check_array, check_random_state
from ..utils import shuffle as util_shuffle
from ..utils.random import sample_without_replacement
from ..externals import six
map = six.moves.map
zip = six.moves.zip


def _generate_hypercube(samples, dimensions, rng):
    """Returns distinct binary samples of length dimensions
    """
    if dimensions > 30:
        return np.hstack([rng.randint(2, size=(samples, dimensions - 30)),
                          _generate_hypercube(samples, 30, rng)])
    out = sample_without_replacement(2 ** dimensions, samples,
                                     random_state=rng).astype(dtype='>u4',
                                                              copy=False)
    out = np.unpackbits(out.view('>u1')).reshape((-1, 32))[:, -dimensions:]
    return out


def make_classification(n_samples=100, n_features=20, n_informative=2,
                        n_redundant=2, n_repeated=0, n_classes=2,
                        n_clusters_per_class=2, weights=None, flip_y=0.01,
                        class_sep=1.0, hypercube=True, shift=0.0, scale=1.0,
                        shuffle=True, random_state=None):
    """Generate a random n-class classification problem.

    This initially creates clusters of points normally distributed (std=1)
    about vertices of an `n_informative`-dimensional hypercube with sides of
    length `2*class_sep` and assigns an equal number of clusters to each
    class. It introduces interdependence between these features and adds
    various types of further noise to the data.

    Prior to shuffling, `X` stacks a number of these primary "informative"
    features, "redundant" linear combinations of these, "repeated" duplicates
    of sampled features, and arbitrary noise for and remaining features.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    n_features : int, optional (default=20)
        The total number of features. These comprise `n_informative`
        informative features, `n_redundant` redundant features, `n_repeated`
        duplicated features and `n_features-n_informative-n_redundant-
        n_repeated` useless features drawn at random.

    n_informative : int, optional (default=2)
        The number of informative features. Each class is composed of a number
        of gaussian clusters each located around the vertices of a hypercube
        in a subspace of dimension `n_informative`. For each cluster,
        informative features are drawn independently from  N(0, 1) and then
        randomly linearly combined within each cluster in order to add
        covariance. The clusters are then placed on the vertices of the
        hypercube.

    n_redundant : int, optional (default=2)
        The number of redundant features. These features are generated as
        random linear combinations of the informative features.

    n_repeated : int, optional (default=0)
        The number of duplicated features, drawn randomly from the informative
        and the redundant features.

    n_classes : int, optional (default=2)
        The number of classes (or labels) of the classification problem.

    n_clusters_per_class : int, optional (default=2)
        The number of clusters per class.

    weights : list of floats or None (default=None)
        The proportions of samples assigned to each class. If None, then
        classes are balanced. Note that if `len(weights) == n_classes - 1`,
        then the last class weight is automatically inferred.
        More than `n_samples` samples may be returned if the sum of `weights`
        exceeds 1.

    flip_y : float, optional (default=0.01)
        The fraction of samples whose class are randomly exchanged. Larger
        values introduce noise in the labels and make the classification
        task harder.

    class_sep : float, optional (default=1.0)
        The factor multiplying the hypercube size.  Larger values spread
        out the clusters/classes and make the classification task easier.

    hypercube : boolean, optional (default=True)
        If True, the clusters are put on the vertices of a hypercube. If
        False, the clusters are put on the vertices of a random polytope.

    shift : float, array of shape [n_features] or None, optional (default=0.0)
        Shift features by the specified value. If None, then features
        are shifted by a random value drawn in [-class_sep, class_sep].

    scale : float, array of shape [n_features] or None, optional (default=1.0)
        Multiply features by the specified value. If None, then features
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

    References
    ----------
    .. [1] I. Guyon, "Design of experiments for the NIPS 2003 variable
           selection benchmark", 2003.

    See also
    --------
    make_blobs: simplified variant
    make_multilabel_classification: unrelated generator for multilabel tasks
    """
    generator = check_random_state(random_state)

    # Count features, clusters and samples
    if n_informative + n_redundant + n_repeated > n_features:
        raise ValueError("Number of informative, redundant and repeated "
                         "features must sum to less than the number of total"
                         " features")
    if 2 ** n_informative < n_classes * n_clusters_per_class:
        raise ValueError("n_classes * n_clusters_per_class must"
                         " be smaller or equal 2 ** n_informative")
    if weights and len(weights) not in [n_classes, n_classes - 1]:
        raise ValueError("Weights specified but incompatible with number "
                         "of classes.")

    n_useless = n_features - n_informative - n_redundant - n_repeated
    n_clusters = n_classes * n_clusters_per_class

    if weights and len(weights) == (n_classes - 1):
        weights = weights + [1.0 - sum(weights)]

    if weights is None:
        weights = [1.0 / n_classes] * n_classes
        weights[-1] = 1.0 - sum(weights[:-1])

    # Distribute samples among clusters by weight
    n_samples_per_cluster = []
    for k in range(n_clusters):
        n_samples_per_cluster.append(int(n_samples * weights[k % n_classes]
                                     / n_clusters_per_class))
    for i in range(n_samples - sum(n_samples_per_cluster)):
        n_samples_per_cluster[i % n_clusters] += 1

    # Initialize X and y
    X = np.zeros((n_samples, n_features))
    y = np.zeros(n_samples, dtype=np.int)

    # Build the polytope whose vertices become cluster centroids
    centroids = _generate_hypercube(n_clusters, n_informative,
                                    generator).astype(float)
    centroids *= 2 * class_sep
    centroids -= class_sep
    if not hypercube:
        centroids *= generator.rand(n_clusters, 1)
        centroids *= generator.rand(1, n_informative)

    # Initially draw informative features from the standard normal
    X[:, :n_informative] = generator.randn(n_samples, n_informative)

    # Create each cluster; a variant of make_blobs
    stop = 0
    for k, centroid in enumerate(centroids):
        start, stop = stop, stop + n_samples_per_cluster[k]
        y[start:stop] = k % n_classes  # assign labels
        X_k = X[start:stop, :n_informative]  # slice a view of the cluster

        A = 2 * generator.rand(n_informative, n_informative) - 1
        X_k[...] = np.dot(X_k, A)  # introduce random covariance

        X_k += centroid  # shift the cluster to a vertex

    # Create redundant features
    if n_redundant > 0:
        B = 2 * generator.rand(n_informative, n_redundant) - 1
        X[:, n_informative:n_informative + n_redundant] = \
            np.dot(X[:, :n_informative], B)

    # Repeat some features
    if n_repeated > 0:
        n = n_informative + n_redundant
        indices = ((n - 1) * generator.rand(n_repeated) + 0.5).astype(np.intp)
        X[:, n:n + n_repeated] = X[:, indices]

    # Fill useless features
    if n_useless > 0:
        X[:, -n_useless:] = generator.randn(n_samples, n_useless)

    # Randomly replace labels
    if flip_y >= 0.0:
        flip_mask = generator.rand(n_samples) < flip_y
        y[flip_mask] = generator.randint(n_classes, size=flip_mask.sum())

    # Randomly shift and scale
    if shift is None:
        shift = (2 * generator.rand(n_features) - 1) * class_sep
    X += shift

    if scale is None:
        scale = 1 + 100 * generator.rand(n_features)
    X *= scale

    if shuffle:
        # Randomly permute samples
        X, y = util_shuffle(X, y, random_state=generator)

        # Randomly permute features
        indices = np.arange(n_features)
        generator.shuffle(indices)
        X[:, :] = X[:, indices]

    return X, y


def make_multilabel_classification(n_samples=100, n_features=20, n_classes=5,
                                   n_labels=2, length=50, allow_unlabeled=True,
                                   sparse=False, return_indicator='dense',
                                   return_distributions=False,
                                   random_state=None):
    """Generate a random multilabel classification problem.

    For each sample, the generative process is:
        - pick the number of labels: n ~ Poisson(n_labels)
        - n times, choose a class c: c ~ Multinomial(theta)
        - pick the document length: k ~ Poisson(length)
        - k times, choose a word: w ~ Multinomial(theta_c)

    In the above process, rejection sampling is used to make sure that
    n is never zero or more than `n_classes`, and that the document length
    is never zero. Likewise, we reject classes which have already been chosen.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    n_features : int, optional (default=20)
        The total number of features.

    n_classes : int, optional (default=5)
        The number of classes of the classification problem.

    n_labels : int, optional (default=2)
        The average number of labels per instance. More precisely, the number
        of labels per sample is drawn from a Poisson distribution with
        ``n_labels`` as its expected value, but samples are bounded (using
        rejection sampling) by ``n_classes``, and must be nonzero if
        ``allow_unlabeled`` is False.

    length : int, optional (default=50)
        The sum of the features (number of words if documents) is drawn from
        a Poisson distribution with this expected value.

    allow_unlabeled : bool, optional (default=True)
        If ``True``, some instances might not belong to any class.

    sparse : bool, optional (default=False)
        If ``True``, return a sparse feature matrix

        .. versionadded:: 0.17
           parameter to allow *sparse* output.

    return_indicator : 'dense' (default) | 'sparse' | False
        If ``dense`` return ``Y`` in the dense binary indicator format. If
        ``'sparse'`` return ``Y`` in the sparse binary indicator format.
        ``False`` returns a list of lists of labels.

    return_distributions : bool, optional (default=False)
        If ``True``, return the prior class probability and conditional
        probabilities of features given classes, from which the data was
        drawn.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, n_features]
        The generated samples.

    Y : array or sparse CSR matrix of shape [n_samples, n_classes]
        The label sets.

    p_c : array, shape [n_classes]
        The probability of each class being drawn. Only returned if
        ``return_distributions=True``.

    p_w_c : array, shape [n_features, n_classes]
        The probability of each feature being drawn given each class.
        Only returned if ``return_distributions=True``.

    """
    generator = check_random_state(random_state)
    p_c = generator.rand(n_classes)
    p_c /= p_c.sum()
    cumulative_p_c = np.cumsum(p_c)
    p_w_c = generator.rand(n_features, n_classes)
    p_w_c /= np.sum(p_w_c, axis=0)

    def sample_example():
        _, n_classes = p_w_c.shape

        # pick a nonzero number of labels per document by rejection sampling
        y_size = n_classes + 1
        while (not allow_unlabeled and y_size == 0) or y_size > n_classes:
            y_size = generator.poisson(n_labels)

        # pick n classes
        y = set()
        while len(y) != y_size:
            # pick a class with probability P(c)
            c = np.searchsorted(cumulative_p_c,
                                generator.rand(y_size - len(y)))
            y.update(c)
        y = list(y)

        # pick a non-zero document length by rejection sampling
        n_words = 0
        while n_words == 0:
            n_words = generator.poisson(length)

        # generate a document of length n_words
        if len(y) == 0:
            # if sample does not belong to any class, generate noise word
            words = generator.randint(n_features, size=n_words)
            return words, y

        # sample words with replacement from selected classes
        cumulative_p_w_sample = p_w_c.take(y, axis=1).sum(axis=1).cumsum()
        cumulative_p_w_sample /= cumulative_p_w_sample[-1]
        words = np.searchsorted(cumulative_p_w_sample, generator.rand(n_words))
        return words, y

    X_indices = array.array('i')
    X_indptr = array.array('i', [0])
    Y = []
    for i in range(n_samples):
        words, y = sample_example()
        X_indices.extend(words)
        X_indptr.append(len(X_indices))
        Y.append(y)
    X_data = np.ones(len(X_indices), dtype=np.float64)
    X = sp.csr_matrix((X_data, X_indices, X_indptr),
                      shape=(n_samples, n_features))
    X.sum_duplicates()
    if not sparse:
        X = X.toarray()

    # return_indicator can be True due to backward compatibility
    if return_indicator in (True, 'sparse', 'dense'):
        lb = MultiLabelBinarizer(sparse_output=(return_indicator == 'sparse'))
        Y = lb.fit([range(n_classes)]).transform(Y)
    elif return_indicator is not False:
        raise ValueError("return_indicator must be either 'sparse', 'dense' "
                         'or False.')
    if return_distributions:
        return X, Y, p_c, p_w_c
    return X, Y


def make_hastie_10_2(n_samples=12000, random_state=None):
    """Generates data for binary classification used in
    Hastie et al. 2009, Example 10.2.

    The ten features are standard independent Gaussian and
    the target ``y`` is defined by::

      y[i] = 1 if np.sum(X[i] ** 2) > 9.34 else -1

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_samples : int, optional (default=12000)
        The number of samples.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, 10]
        The input samples.

    y : array of shape [n_samples]
        The output values.

    References
    ----------
    .. [1] T. Hastie, R. Tibshirani and J. Friedman, "Elements of Statistical
           Learning Ed. 2", Springer, 2009.

    See also
    --------
    make_gaussian_quantiles: a generalization of this dataset approach
    """
    rs = check_random_state(random_state)

    shape = (n_samples, 10)
    X = rs.normal(size=shape).reshape(shape)
    y = ((X ** 2.0).sum(axis=1) > 9.34).astype(np.float64)
    y[y == 0.0] = -1.0

    return X, y


def make_regression(n_samples=100, n_features=100, n_informative=10,
                    n_targets=1, bias=0.0, effective_rank=None,
                    tail_strength=0.5, noise=0.0, shuffle=True, coef=False,
                    random_state=None):
    """Generate a random regression problem.

    The input set can either be well conditioned (by default) or have a low
    rank-fat tail singular profile. See :func:`make_low_rank_matrix` for
    more details.

    The output is generated by applying a (potentially biased) random linear
    regression model with `n_informative` nonzero regressors to the previously
    generated input and some gaussian centered noise with some adjustable
    scale.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The number of samples.

    n_features : int, optional (default=100)
        The number of features.

    n_informative : int, optional (default=10)
        The number of informative features, i.e., the number of features used
        to build the linear model used to generate the output.

    n_targets : int, optional (default=1)
        The number of regression targets, i.e., the dimension of the y output
        vector associated with a sample. By default, the output is a scalar.

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

    y : array of shape [n_samples] or [n_samples, n_targets]
        The output values.

    coef : array of shape [n_features] or [n_features, n_targets], optional
        The coefficient of the underlying linear model. It is returned only if
        coef is True.
    """
    n_informative = min(n_features, n_informative)
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
    ground_truth = np.zeros((n_features, n_targets))
    ground_truth[:n_informative, :] = 100 * generator.rand(n_informative,
                                                           n_targets)

    y = np.dot(X, ground_truth) + bias

    # Add noise
    if noise > 0.0:
        y += generator.normal(scale=noise, size=y.shape)

    # Randomly permute samples and features
    if shuffle:
        X, y = util_shuffle(X, y, random_state=generator)

        indices = np.arange(n_features)
        generator.shuffle(indices)
        X[:, :] = X[:, indices]
        ground_truth = ground_truth[indices]

    y = np.squeeze(y)

    if coef:
        return X, y, np.squeeze(ground_truth)

    else:
        return X, y


def make_circles(n_samples=100, shuffle=True, noise=None, random_state=None,
                 factor=.8):
    """Make a large circle containing a smaller circle in 2d.

    A simple toy dataset to visualize clustering and classification
    algorithms.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The total number of points generated.

    shuffle : bool, optional (default=True)
        Whether to shuffle the samples.

    noise : double or None (default=None)
        Standard deviation of Gaussian noise added to the data.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    factor : double < 1 (default=.8)
        Scale factor between inner and outer circle.

    Returns
    -------
    X : array of shape [n_samples, 2]
        The generated samples.

    y : array of shape [n_samples]
        The integer labels (0 or 1) for class membership of each sample.
    """

    if factor > 1 or factor < 0:
        raise ValueError("'factor' has to be between 0 and 1.")

    generator = check_random_state(random_state)
    # so as not to have the first point = last point, we add one and then
    # remove it.
    linspace = np.linspace(0, 2 * np.pi, n_samples // 2 + 1)[:-1]
    outer_circ_x = np.cos(linspace)
    outer_circ_y = np.sin(linspace)
    inner_circ_x = outer_circ_x * factor
    inner_circ_y = outer_circ_y * factor

    X = np.vstack((np.append(outer_circ_x, inner_circ_x),
                   np.append(outer_circ_y, inner_circ_y))).T
    y = np.hstack([np.zeros(n_samples // 2, dtype=np.intp),
                   np.ones(n_samples // 2, dtype=np.intp)])
    if shuffle:
        X, y = util_shuffle(X, y, random_state=generator)

    if noise is not None:
        X += generator.normal(scale=noise, size=X.shape)

    return X, y


def make_moons(n_samples=100, shuffle=True, noise=None, random_state=None):
    """Make two interleaving half circles

    A simple toy dataset to visualize clustering and classification
    algorithms. Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The total number of points generated.

    shuffle : bool, optional (default=True)
        Whether to shuffle the samples.

    noise : double or None (default=None)
        Standard deviation of Gaussian noise added to the data.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape [n_samples, 2]
        The generated samples.

    y : array of shape [n_samples]
        The integer labels (0 or 1) for class membership of each sample.
    """

    n_samples_out = n_samples // 2
    n_samples_in = n_samples - n_samples_out

    generator = check_random_state(random_state)

    outer_circ_x = np.cos(np.linspace(0, np.pi, n_samples_out))
    outer_circ_y = np.sin(np.linspace(0, np.pi, n_samples_out))
    inner_circ_x = 1 - np.cos(np.linspace(0, np.pi, n_samples_in))
    inner_circ_y = 1 - np.sin(np.linspace(0, np.pi, n_samples_in)) - .5

    X = np.vstack((np.append(outer_circ_x, inner_circ_x),
                   np.append(outer_circ_y, inner_circ_y))).T
    y = np.hstack([np.zeros(n_samples_out, dtype=np.intp),
                   np.ones(n_samples_in, dtype=np.intp)])

    if shuffle:
        X, y = util_shuffle(X, y, random_state=generator)

    if noise is not None:
        X += generator.normal(scale=noise, size=X.shape)

    return X, y


def make_blobs(n_samples=100, n_features=2, centers=3, cluster_std=1.0,
               center_box=(-10.0, 10.0), shuffle=True, random_state=None):
    """Generate isotropic Gaussian blobs for clustering.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_samples : int, optional (default=100)
        The total number of points equally divided among clusters.

    n_features : int, optional (default=2)
        The number of features for each sample.

    centers : int or array of shape [n_centers, n_features], optional
        (default=3)
        The number of centers to generate, or the fixed center locations.

    cluster_std : float or sequence of floats, optional (default=1.0)
        The standard deviation of the clusters.

    center_box : pair of floats (min, max), optional (default=(-10.0, 10.0))
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
    >>> print(X.shape)
    (10, 2)
    >>> y
    array([0, 0, 1, 0, 2, 2, 2, 1, 1, 0])

    See also
    --------
    make_classification: a more intricate variant
    """
    generator = check_random_state(random_state)

    if isinstance(centers, numbers.Integral):
        centers = generator.uniform(center_box[0], center_box[1],
                                    size=(centers, n_features))
    else:
        centers = check_array(centers)
        n_features = centers.shape[1]

    if isinstance(cluster_std, numbers.Real):
        cluster_std = np.ones(len(centers)) * cluster_std

    X = []
    y = []

    n_centers = centers.shape[0]
    n_samples_per_center = [int(n_samples // n_centers)] * n_centers

    for i in range(n_samples % n_centers):
        n_samples_per_center[i] += 1

    for i, (n, std) in enumerate(zip(n_samples_per_center, cluster_std)):
        X.append(centers[i] + generator.normal(scale=std,
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
    """Generate the "Friedman \#1" regression problem

    This dataset is described in Friedman [1] and Breiman [2].

    Inputs `X` are independent features uniformly distributed on the interval
    [0, 1]. The output `y` is created according to the formula::

        y(X) = 10 * sin(pi * X[:, 0] * X[:, 1]) + 20 * (X[:, 2] - 0.5) ** 2 \
+ 10 * X[:, 3] + 5 * X[:, 4] + noise * N(0, 1).

    Out of the `n_features` features, only 5 are actually used to compute
    `y`. The remaining features are independent of `y`.

    The number of features has to be >= 5.

    Read more in the :ref:`User Guide <sample_generators>`.

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

    References
    ----------
    .. [1] J. Friedman, "Multivariate adaptive regression splines", The Annals
           of Statistics 19 (1), pages 1-67, 1991.

    .. [2] L. Breiman, "Bagging predictors", Machine Learning 24,
           pages 123-140, 1996.
    """
    if n_features < 5:
        raise ValueError("n_features must be at least five.")

    generator = check_random_state(random_state)

    X = generator.rand(n_samples, n_features)
    y = 10 * np.sin(np.pi * X[:, 0] * X[:, 1]) + 20 * (X[:, 2] - 0.5) ** 2 \
        + 10 * X[:, 3] + 5 * X[:, 4] + noise * generator.randn(n_samples)

    return X, y


def make_friedman2(n_samples=100, noise=0.0, random_state=None):
    """Generate the "Friedman \#2" regression problem

    This dataset is described in Friedman [1] and Breiman [2].

    Inputs `X` are 4 independent features uniformly distributed on the
    intervals::

        0 <= X[:, 0] <= 100,
        40 * pi <= X[:, 1] <= 560 * pi,
        0 <= X[:, 2] <= 1,
        1 <= X[:, 3] <= 11.

    The output `y` is created according to the formula::

        y(X) = (X[:, 0] ** 2 + (X[:, 1] * X[:, 2] \
 - 1 / (X[:, 1] * X[:, 3])) ** 2) ** 0.5 + noise * N(0, 1).

    Read more in the :ref:`User Guide <sample_generators>`.

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

    References
    ----------
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
    """Generate the "Friedman \#3" regression problem

    This dataset is described in Friedman [1] and Breiman [2].

    Inputs `X` are 4 independent features uniformly distributed on the
    intervals::

        0 <= X[:, 0] <= 100,
        40 * pi <= X[:, 1] <= 560 * pi,
        0 <= X[:, 2] <= 1,
        1 <= X[:, 3] <= 11.

    The output `y` is created according to the formula::

        y(X) = arctan((X[:, 1] * X[:, 2] - 1 / (X[:, 1] * X[:, 3])) \
/ X[:, 0]) + noise * N(0, 1).

    Read more in the :ref:`User Guide <sample_generators>`.

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

    References
    ----------
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

    Read more in the :ref:`User Guide <sample_generators>`.

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
    u, _ = linalg.qr(generator.randn(n_samples, n), mode='economic')
    v, _ = linalg.qr(generator.randn(n_features, n), mode='economic')

    # Index of the singular values
    singular_ind = np.arange(n, dtype=np.float64)

    # Build the singular profile by assembling signal and noise components
    low_rank = ((1 - tail_strength) *
                np.exp(-1.0 * (singular_ind / effective_rank) ** 2))
    tail = tail_strength * np.exp(-0.1 * singular_ind / effective_rank)
    s = np.identity(n) * (low_rank + tail)

    return np.dot(np.dot(u, s), v.T)


def make_sparse_coded_signal(n_samples, n_components, n_features,
                             n_nonzero_coefs, random_state=None):
    """Generate a signal as a sparse combination of dictionary elements.

    Returns a matrix Y = DX, such as D is (n_features, n_components),
    X is (n_components, n_samples) and each column of X has exactly
    n_nonzero_coefs non-zero elements.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    n_samples : int
        number of samples to generate

    n_components :  int,
        number of components in the dictionary

    n_features : int
        number of features of the dataset to generate

    n_nonzero_coefs : int
        number of active (non-zero) coefficients in each sample

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    data : array of shape [n_features, n_samples]
        The encoded signal (Y).

    dictionary : array of shape [n_features, n_components]
        The dictionary with normalized components (D).

    code : array of shape [n_components, n_samples]
        The sparse code such that each column of this matrix has exactly
        n_nonzero_coefs non-zero items (X).

    """
    generator = check_random_state(random_state)

    # generate dictionary
    D = generator.randn(n_features, n_components)
    D /= np.sqrt(np.sum((D ** 2), axis=0))

    # generate code
    X = np.zeros((n_components, n_samples))
    for i in range(n_samples):
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

    Read more in the :ref:`User Guide <sample_generators>`.

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

    References
    ----------
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

    Read more in the :ref:`User Guide <sample_generators>`.

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

    See also
    --------
    make_sparse_spd_matrix
    """
    generator = check_random_state(random_state)

    A = generator.rand(n_dim, n_dim)
    U, s, V = linalg.svd(np.dot(A.T, A))
    X = np.dot(np.dot(U, 1.0 + np.diag(generator.rand(n_dim))), V)

    return X


def make_sparse_spd_matrix(dim=1, alpha=0.95, norm_diag=False,
                           smallest_coef=.1, largest_coef=.9,
                           random_state=None):
    """Generate a sparse symmetric definite positive matrix.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    dim : integer, optional (default=1)
        The size of the random matrix to generate.

    alpha : float between 0 and 1, optional (default=0.95)
        The probability that a coefficient is zero (see notes). Larger values
        enforce more sparsity.

    norm_diag : boolean, optional (default=False)
        Whether to normalize the output matrix to make the leading diagonal
        elements all 1

    smallest_coef : float between 0 and 1, optional (default=0.1)
        The value of the smallest coefficient.

    largest_coef : float between 0 and 1, optional (default=0.9)
        The value of the largest coefficient.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    prec : sparse matrix of shape (dim, dim)
        The generated matrix.

    Notes
    -----
    The sparsity is actually imposed on the cholesky factor of the matrix.
    Thus alpha does not translate directly into the filling fraction of
    the matrix itself.

    See also
    --------
    make_spd_matrix
    """
    random_state = check_random_state(random_state)

    chol = -np.eye(dim)
    aux = random_state.rand(dim, dim)
    aux[aux < alpha] = 0
    aux[aux > alpha] = (smallest_coef
                        + (largest_coef - smallest_coef)
                        * random_state.rand(np.sum(aux > alpha)))
    aux = np.tril(aux, k=-1)

    # Permute the lines: we don't want to have asymmetries in the final
    # SPD matrix
    permutation = random_state.permutation(dim)
    aux = aux[permutation].T[permutation]
    chol += aux
    prec = np.dot(chol.T, chol)

    if norm_diag:
        # Form the diagonal vector into a row matrix
        d = np.diag(prec).reshape(1, prec.shape[0])
        d = 1. / np.sqrt(d)

        prec *= d
        prec *= d.T

    return prec


def make_swiss_roll(n_samples=100, noise=0.0, random_state=None):
    """Generate a swiss roll dataset.

    Read more in the :ref:`User Guide <sample_generators>`.

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

    References
    ----------
    .. [1] S. Marsland, "Machine Learning: An Algorithmic Perspective",
           Chapter 10, 2009.
           http://seat.massey.ac.nz/personal/s.r.marsland/Code/10/lle.py
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

    Read more in the :ref:`User Guide <sample_generators>`.

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


def make_gaussian_quantiles(mean=None, cov=1., n_samples=100,
                            n_features=2, n_classes=3,
                            shuffle=True, random_state=None):
    """Generate isotropic Gaussian and label samples by quantile

    This classification dataset is constructed by taking a multi-dimensional
    standard normal distribution and defining classes separated by nested
    concentric multi-dimensional spheres such that roughly equal numbers of
    samples are in each class (quantiles of the :math:`\chi^2` distribution).

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    mean : array of shape [n_features], optional (default=None)
        The mean of the multi-dimensional normal distribution.
        If None then use the origin (0, 0, ...).

    cov : float, optional (default=1.)
        The covariance matrix will be this value times the unit matrix. This
        dataset only produces symmetric normal distributions.

    n_samples : int, optional (default=100)
        The total number of points equally divided among classes.

    n_features : int, optional (default=2)
        The number of features for each sample.

    n_classes : int, optional (default=3)
        The number of classes

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
        The integer labels for quantile membership of each sample.

    Notes
    -----
    The dataset is from Zhu et al [1].

    References
    ----------
    .. [1] J. Zhu, H. Zou, S. Rosset, T. Hastie, "Multi-class AdaBoost", 2009.

    """
    if n_samples < n_classes:
        raise ValueError("n_samples must be at least n_classes")

    generator = check_random_state(random_state)

    if mean is None:
        mean = np.zeros(n_features)
    else:
        mean = np.array(mean)

    # Build multivariate normal distribution
    X = generator.multivariate_normal(mean, cov * np.identity(n_features),
                                      (n_samples,))

    # Sort by distance from origin
    idx = np.argsort(np.sum((X - mean[np.newaxis, :]) ** 2, axis=1))
    X = X[idx, :]

    # Label by quantile
    step = n_samples // n_classes

    y = np.hstack([np.repeat(np.arange(n_classes), step),
                   np.repeat(n_classes - 1, n_samples - step * n_classes)])

    if shuffle:
        X, y = util_shuffle(X, y, random_state=generator)

    return X, y


def _shuffle(data, random_state=None):
    generator = check_random_state(random_state)
    n_rows, n_cols = data.shape
    row_idx = generator.permutation(n_rows)
    col_idx = generator.permutation(n_cols)
    result = data[row_idx][:, col_idx]
    return result, row_idx, col_idx


def make_biclusters(shape, n_clusters, noise=0.0, minval=10,
                    maxval=100, shuffle=True, random_state=None):
    """Generate an array with constant block diagonal structure for
    biclustering.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    shape : iterable (n_rows, n_cols)
        The shape of the result.

    n_clusters : integer
        The number of biclusters.

    noise : float, optional (default=0.0)
        The standard deviation of the gaussian noise.

    minval : int, optional (default=10)
        Minimum value of a bicluster.

    maxval : int, optional (default=100)
        Maximum value of a bicluster.

    shuffle : boolean, optional (default=True)
        Shuffle the samples.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape `shape`
        The generated array.

    rows : array of shape (n_clusters, X.shape[0],)
        The indicators for cluster membership of each row.

    cols : array of shape (n_clusters, X.shape[1],)
        The indicators for cluster membership of each column.

    References
    ----------

    .. [1] Dhillon, I. S. (2001, August). Co-clustering documents and
        words using bipartite spectral graph partitioning. In Proceedings
        of the seventh ACM SIGKDD international conference on Knowledge
        discovery and data mining (pp. 269-274). ACM.

    See also
    --------
    make_checkerboard
    """
    generator = check_random_state(random_state)
    n_rows, n_cols = shape
    consts = generator.uniform(minval, maxval, n_clusters)

    # row and column clusters of approximately equal sizes
    row_sizes = generator.multinomial(n_rows,
                                      np.repeat(1.0 / n_clusters,
                                                n_clusters))
    col_sizes = generator.multinomial(n_cols,
                                      np.repeat(1.0 / n_clusters,
                                                n_clusters))

    row_labels = np.hstack(list(np.repeat(val, rep) for val, rep in
                                zip(range(n_clusters), row_sizes)))
    col_labels = np.hstack(list(np.repeat(val, rep) for val, rep in
                                zip(range(n_clusters), col_sizes)))

    result = np.zeros(shape, dtype=np.float64)
    for i in range(n_clusters):
        selector = np.outer(row_labels == i, col_labels == i)
        result[selector] += consts[i]

    if noise > 0:
        result += generator.normal(scale=noise, size=result.shape)

    if shuffle:
        result, row_idx, col_idx = _shuffle(result, random_state)
        row_labels = row_labels[row_idx]
        col_labels = col_labels[col_idx]

    rows = np.vstack(row_labels == c for c in range(n_clusters))
    cols = np.vstack(col_labels == c for c in range(n_clusters))

    return result, rows, cols


def make_checkerboard(shape, n_clusters, noise=0.0, minval=10,
                      maxval=100, shuffle=True, random_state=None):

    """Generate an array with block checkerboard structure for
    biclustering.

    Read more in the :ref:`User Guide <sample_generators>`.

    Parameters
    ----------
    shape : iterable (n_rows, n_cols)
        The shape of the result.

    n_clusters : integer or iterable (n_row_clusters, n_column_clusters)
        The number of row and column clusters.

    noise : float, optional (default=0.0)
        The standard deviation of the gaussian noise.

    minval : int, optional (default=10)
        Minimum value of a bicluster.

    maxval : int, optional (default=100)
        Maximum value of a bicluster.

    shuffle : boolean, optional (default=True)
        Shuffle the samples.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Returns
    -------
    X : array of shape `shape`
        The generated array.

    rows : array of shape (n_clusters, X.shape[0],)
        The indicators for cluster membership of each row.

    cols : array of shape (n_clusters, X.shape[1],)
        The indicators for cluster membership of each column.


    References
    ----------

    .. [1] Kluger, Y., Basri, R., Chang, J. T., & Gerstein, M. (2003).
        Spectral biclustering of microarray data: coclustering genes
        and conditions. Genome research, 13(4), 703-716.

    See also
    --------
    make_biclusters
    """
    generator = check_random_state(random_state)

    if hasattr(n_clusters, "__len__"):
        n_row_clusters, n_col_clusters = n_clusters
    else:
        n_row_clusters = n_col_clusters = n_clusters

    # row and column clusters of approximately equal sizes
    n_rows, n_cols = shape
    row_sizes = generator.multinomial(n_rows,
                                      np.repeat(1.0 / n_row_clusters,
                                                n_row_clusters))
    col_sizes = generator.multinomial(n_cols,
                                      np.repeat(1.0 / n_col_clusters,
                                                n_col_clusters))

    row_labels = np.hstack(list(np.repeat(val, rep) for val, rep in
                                zip(range(n_row_clusters), row_sizes)))
    col_labels = np.hstack(list(np.repeat(val, rep) for val, rep in
                                zip(range(n_col_clusters), col_sizes)))

    result = np.zeros(shape, dtype=np.float64)
    for i in range(n_row_clusters):
        for j in range(n_col_clusters):
            selector = np.outer(row_labels == i, col_labels == j)
            result[selector] += generator.uniform(minval, maxval)

    if noise > 0:
        result += generator.normal(scale=noise, size=result.shape)

    if shuffle:
        result, row_idx, col_idx = _shuffle(result, random_state)
        row_labels = row_labels[row_idx]
        col_labels = col_labels[col_idx]

    rows = np.vstack(row_labels == label
                     for label in range(n_row_clusters)
                     for _ in range(n_col_clusters))
    cols = np.vstack(col_labels == label
                     for _ in range(n_row_clusters)
                     for label in range(n_col_clusters))

    return result, rows, cols
