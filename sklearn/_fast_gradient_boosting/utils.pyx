# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
"""This module contains utility routines."""

from cython.parallel import prange

from .binning import BinMapper
from .types cimport Y_DTYPE_C


def get_lightgbm_estimator(pygbm_estimator):
    """Return an unfitted LightGBM estimator with matching hyperparams.

    This utility function takes care of renaming the PyGBM parameters into
    their LightGBM equivalent parameters.
    """
    from lightgbm import LGBMRegressor
    from lightgbm import LGBMClassifier

    # Import here to avoid cyclic dependencies
    from .gradient_boosting import FastGradientBoostingClassifier

    pygbm_params = pygbm_estimator.get_params()

    if pygbm_params['loss'] == 'auto':
        raise ValueError('auto loss is not accepted. We need to know if '
                         'the problem is binary or multiclass classification.')
    if pygbm_params['n_iter_no_change'] is not None:
        raise NotImplementedError('Early stopping should be deactivated.')

    loss_mapping = {
        'least_squares': 'regression_l2',
        'binary_crossentropy': 'binary',
        'categorical_crossentropy': 'multiclass'
    }

    lgbm_params = {
        'objective': loss_mapping[pygbm_params['loss']],
        'learning_rate': pygbm_params['learning_rate'],
        'n_estimators': pygbm_params['n_estimators'],
        'num_leaves': pygbm_params['max_leaf_nodes'],
        'max_depth': pygbm_params['max_depth'],
        'min_child_samples': pygbm_params['min_samples_leaf'],
        'reg_lambda': pygbm_params['l2_regularization'],
        'max_bin': pygbm_params['max_bins'],
        'min_data_in_bin': 1,
        'min_child_weight': 1e-3,
        'min_sum_hessian_in_leaf': 1e-3,
        'min_split_gain': 0,
        'verbosity': 10 if pygbm_params['verbose'] else -10,
        'boost_from_average': True,
        'enable_bundle': False,  # also makes feature order consistent
        'min_data_in_bin': 1,
        'subsample_for_bin': BinMapper().subsample,
    }
    # TODO: change hardcoded values when / if they're arguments to the
    # estimator.

    if pygbm_params['loss'] == 'categorical_crossentropy':
        # LGBM multiplies hessians by 2 in multiclass loss.
        lgbm_params['min_sum_hessian_in_leaf'] *= 2
        lgbm_params['learning_rate'] *= 2

    if isinstance(pygbm_estimator, FastGradientBoostingClassifier):
        Est = LGBMClassifier
    else:
        Est = LGBMRegressor

    return Est(**lgbm_params)


def sum_parallel(Y_DTYPE_C [:] array):

    cdef:
        Y_DTYPE_C out = 0.
        int i = 0

    with nogil:
        for i in prange(array.shape[0], schedule='static'):
            out += array[i]

    return out
