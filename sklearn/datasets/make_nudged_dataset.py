"""
Generate images that are nudged to augment images dataset
"""

# Authors: G. Synnaeve
# License: BSD 3 clause

from sklearn import datasets
# TODO example/test with datasets.load_digits() (.data np.asarray / .target)

import numpy as np
from scipy.ndimage import convolve

from ..utils import check_random_state


def nudge_dataset(X, Y, images_shape=(8, 8, 1),
        n_samples=100, nudging_strength=1, random_state=None):
    """ 
    Generate randomly nudged (left-right-up-down) images of the the ones in
    the given dataset.

    X : array of shape [n_examples, n_features]
        images examples

    y : array of shape [n_examples] 
        images labels

    images_shape : 3 elements tuple, (images_x_size, images_y_size, #colors)
        (default=(8,8,1))

    n_samples : number of additional samples generated from the dataset
        (default=100)

    nudging_strength : number of pixels by which to nudge the images
        (default=1)

    random_state : random state initialization
        (default=None)

    Returns
    -------

    X : array of shape [n_examples + n_samples, n_features]
        images examples

    y : array of shape [n_examples + n_samples] 
        images labels
    """
    rng = check_random_state(random_state)

    direction_vectors = [
        [[0, nudging_strength, 0],
         [0, 0, 0],
         [0, 0, 0]],

        [[0, 0, 0],
         [nudging_strength, 0, 0],
         [0, 0, 0]],

        [[0, 0, 0],
         [0, 0, nudging_strength],
         [0, 0, 0]],

        [[0, 0, 0],
         [0, 0, 0],
         [0, nudging_strength, 0]]] # TODO make that a parameter

    shift = lambda x, w: convolve(x.reshape((8, 8)), mode='constant',
                                  weights=w).ravel()
    sample_dir = None # TODO
    # TODO scipy.ndimage.rotate(img, angle_in_deg, reshape=False)

    i_ = rng.randint(0, X.shape[0], n_samples)
    subset_X = X[i_]
    X = np.concatenate([X] +
                       [np.apply_along_axis(shift, 1, subset_X, vector)
                        for vector in sample_dir(direction_vectors)])
    y = np.concatenate([y] + [y[i_]], axis=0)
    return X, y

