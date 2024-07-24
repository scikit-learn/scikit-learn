"""This module contains utility routines."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ...base import is_classifier
from .binning import _BinMapper


def get_equivalent_estimator(estimator, lib="lightgbm", n_classes=None):
    """Return an unfitted estimator from another lib with matching hyperparams.

    This utility function takes care of renaming the sklearn parameters into
    their LightGBM, XGBoost or CatBoost equivalent parameters.

    # unmapped XGB parameters:
    # - min_samples_leaf
    # - min_data_in_bin
    # - min_split_gain (there is min_split_loss though?)

    # unmapped Catboost parameters:
    # max_leaves
    # min_*
    """

    if lib not in ("lightgbm", "xgboost", "catboost"):
        raise ValueError(
            "accepted libs are lightgbm, xgboost, and catboost.  got {}".format(lib)
        )

    sklearn_params = estimator.get_params()

    if sklearn_params["loss"] == "auto":
        raise ValueError(
            "auto loss is not accepted. We need to know if "
            "the problem is binary or multiclass classification."
        )
    if sklearn_params["early_stopping"]:
        raise NotImplementedError("Early stopping should be deactivated.")

    lightgbm_loss_mapping = {
        "squared_error": "regression_l2",
        "absolute_error": "regression_l1",
        "log_loss": "binary" if n_classes == 2 else "multiclass",
        "gamma": "gamma",
        "poisson": "poisson",
    }

    lightgbm_params = {
        "objective": lightgbm_loss_mapping[sklearn_params["loss"]],
        "learning_rate": sklearn_params["learning_rate"],
        "n_estimators": sklearn_params["max_iter"],
        "num_leaves": sklearn_params["max_leaf_nodes"],
        "max_depth": sklearn_params["max_depth"],
        "min_data_in_leaf": sklearn_params["min_samples_leaf"],
        "reg_lambda": sklearn_params["l2_regularization"],
        "max_bin": sklearn_params["max_bins"],
        "min_data_in_bin": 1,
        "min_sum_hessian_in_leaf": 1e-3,
        "min_split_gain": 0,
        "verbosity": 10 if sklearn_params["verbose"] else -10,
        "boost_from_average": True,
        "enable_bundle": False,  # also makes feature order consistent
        "subsample_for_bin": _BinMapper().subsample,
        "poisson_max_delta_step": 1e-12,
        "feature_fraction_bynode": sklearn_params["max_features"],
    }

    if sklearn_params["loss"] == "log_loss" and n_classes > 2:
        # LightGBM multiplies hessians by 2 in multiclass loss.
        lightgbm_params["min_sum_hessian_in_leaf"] *= 2
        # LightGBM 3.0 introduced a different scaling of the hessian for the multiclass
        # case.
        # It is equivalent of scaling the learning rate.
        # See https://github.com/microsoft/LightGBM/pull/3256.
        if n_classes is not None:
            lightgbm_params["learning_rate"] *= n_classes / (n_classes - 1)

    # XGB
    xgboost_loss_mapping = {
        "squared_error": "reg:linear",
        "absolute_error": "LEAST_ABSOLUTE_DEV_NOT_SUPPORTED",
        "log_loss": "reg:logistic" if n_classes == 2 else "multi:softmax",
        "gamma": "reg:gamma",
        "poisson": "count:poisson",
    }

    xgboost_params = {
        "tree_method": "hist",
        "grow_policy": "lossguide",  # so that we can set max_leaves
        "objective": xgboost_loss_mapping[sklearn_params["loss"]],
        "learning_rate": sklearn_params["learning_rate"],
        "n_estimators": sklearn_params["max_iter"],
        "max_leaves": sklearn_params["max_leaf_nodes"],
        "max_depth": sklearn_params["max_depth"] or 0,
        "lambda": sklearn_params["l2_regularization"],
        "max_bin": sklearn_params["max_bins"],
        "min_child_weight": 1e-3,
        "verbosity": 2 if sklearn_params["verbose"] else 0,
        "silent": sklearn_params["verbose"] == 0,
        "n_jobs": -1,
        "colsample_bynode": sklearn_params["max_features"],
    }

    # Catboost
    catboost_loss_mapping = {
        "squared_error": "RMSE",
        # catboost does not support MAE when leaf_estimation_method is Newton
        "absolute_error": "LEAST_ASBOLUTE_DEV_NOT_SUPPORTED",
        "log_loss": "Logloss" if n_classes == 2 else "MultiClass",
        "gamma": None,
        "poisson": "Poisson",
    }

    catboost_params = {
        "loss_function": catboost_loss_mapping[sklearn_params["loss"]],
        "learning_rate": sklearn_params["learning_rate"],
        "iterations": sklearn_params["max_iter"],
        "depth": sklearn_params["max_depth"],
        "reg_lambda": sklearn_params["l2_regularization"],
        "max_bin": sklearn_params["max_bins"],
        "feature_border_type": "Median",
        "leaf_estimation_method": "Newton",
        "verbose": bool(sklearn_params["verbose"]),
    }

    if lib == "lightgbm":
        from lightgbm import LGBMClassifier, LGBMRegressor

        if is_classifier(estimator):
            return LGBMClassifier(**lightgbm_params)
        else:
            return LGBMRegressor(**lightgbm_params)

    elif lib == "xgboost":
        from xgboost import XGBClassifier, XGBRegressor

        if is_classifier(estimator):
            return XGBClassifier(**xgboost_params)
        else:
            return XGBRegressor(**xgboost_params)

    else:
        from catboost import CatBoostClassifier, CatBoostRegressor

        if is_classifier(estimator):
            return CatBoostClassifier(**catboost_params)
        else:
            return CatBoostRegressor(**catboost_params)
