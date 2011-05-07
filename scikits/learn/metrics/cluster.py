"""Utilities to evaluate the clustering performance of models

Functions named as *_score return a scalar value to maximize: the higher the
better
"""

# Authors: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD Style.

import numpy as np


def homogeneity_completeness(labels_true, labels_pred):
    """Compute the pair of homogeneity and completeness scores at once

    Those metrics are based on normalized conditional entropy measures of
    the clustering labeling to evaluate given the knowledge of a Ground
    Truth class labels of the same samples.

    A clustering result satisfies homogeneity if all of its clusters
    contain only data points which are members of a single class.

    A clustering result satisfies completeness if all the data points
    that are members of a given class are elements of the same cluster.

    Both scores have positive values between 0.0 and 1.0, larger values
    being desirable.

    Parameters
    ----------
    labels_true : int array, shape = [n_samples]
        ground truth class labels to be used as a reference

    labels_pred : array, shape = [n_samples]
        cluster labels to evaluate

    Returns
    -------
    homogeneity: float
       score between 0.0 and 1.0. 1.0 stands for perfectly homogenous labeling

    completeness: float
       score between 0.0 and 1.0. 1.0 stands for perfectly complete labeling

    See also
    --------
    - homogeneity_score
    - completeness_score
    - v_measure_score
    """
    classes = np.unique(labels_true)
    clusters = np.unique(labels_pred)


def homogeneity_score(labels_true, labels_pred):
    """Homogeneity metric of a cluster labeling given a ground truth

    A clustering result satisfies homogeneity if all of its clusters
    contain only data points which are members of a single class.

    Parameters
    ----------
    labels_true : int array, shape = [n_samples]
        ground truth class labels to be used as a reference

    labels_pred : array, shape = [n_samples]
        cluster labels to evaluate

    Returns
    -------
    homogeneity: float
       score between 0.0 and 1.0. 1.0 stands for perfectly homogenous labeling

    References
    ----------
    V-Measure: A conditional entropy-based external cluster evaluation measure
    Andrew Rosenberg and Julia Hirschberg, 2007
    http://acl.ldc.upenn.edu/D/D07/D07-1043.pdf

    See also
    --------
    - completeness_score
    - v_measure_score

    Examples
    --------

    Perfect labelings are homegenous::

      >>> homogeneity_score([0, 0, 1, 1], [1, 1, 0, 0])
      1.0

    Non-pefect labelings that futher split classes into more clusters can be
    perfectly homogenous:

      >>> homogeneity_score([0, 0, 1, 1], [0, 0, 1, 2])
      1.0
      >>> homogeneity_score([0, 0, 1, 1], [0, 1, 2, 3])
      1.0

    Clusters that include samples from different classes do not homogenous
    labeling::

      >>> homogeneity_score([0, 0, 1, 1], [0, 1, 0, 1])
      0.0
      >>> homogeneity_score([0, 0, 1, 1], [0, 0, 0, 0])
      0.0

    """
    return homogeneity_completeness(labels_true, labels_pred)[0]


def completeness_score(labels_true, labels_pred):
    """Homogeneity metric of a cluster labeling given a ground truth

    A clustering result satisfies completeness if all the data points
    that are members of a given class are elements of the same cluster.

    Parameters
    ----------
    labels_true : int array, shape = [n_samples]
        ground truth class labels to be used as a reference

    labels_pred : array, shape = [n_samples]
        cluster labels to evaluate

    Returns
    -------
    completeness: float
       score between 0.0 and 1.0. 1.0 stands for perfectly complete labeling

    References
    ----------
    V-Measure: A conditional entropy-based external cluster evaluation measure
    Andrew Rosenberg and Julia Hirschberg, 2007
    http://acl.ldc.upenn.edu/D/D07/D07-1043.pdf

    See also
    --------
    - homogeneity_score
    - v_measure_score

    Examples
    --------

    >>> # TODO

    """
    return homogeneity_completeness(labels_true, labels_pred)[1]
