"""
==============================================
Supporting callbacks in third party estimators
==============================================

.. currentmodule:: sklearn

This document shows how to make third party :term:`estimators` and
:term:`meta-estimators` compatible with the callback mechanisms supported by
scikit-learn.

Generally speaking, a callback is a function that is provided by the user to be called
at specific steps of a process, or to be triggered by specific events. Callbacks provide
a clean mechanism for inserting custom logic like monitoring progress or metrics,
without modifying the core algorithm of the process.

In scikit-learn, callbacks take the form of classes following a protocol. This protocol
requires the callback classes to implement specific methods (referred as callback hooks)
which will be called at specific steps of the fitting of an estimator or a
meta-estimator. These specific methods are
:meth:`~callback._base.Callback.on_fit_begin`,
:meth:`~callback._base.Callback.on_fit_task_end` and
:meth:`~callback._base.Callback.on_fit_end`, which are respectively called at the start
of the :term:`fit` method, at the end of each task in ``fit`` and at the end of the
``fit`` method. In scikit-learn estimators, a ``fit`` task is usually one step of a
loop, with nested loops corresponding to netsed tasks. In general a task can be whatever
unit of work the estimator's developer wants it to be, however if the defined tasks are
too unusual it might comprommise the compatibility of the estimator with scikit-learn's
builtin callbacks.

In order to support the callbacks, estimators need to initialize and manage
:class:`~callback.CallbackContext` objects. As the name implies, these objects hold the
contextual information necessary to run the callback hooks. They are also responsible
for calling the callback hooks at the right time.


In the following, we will show how to convert an example estimator class and an example
meta-estimator class to make them compliant with the scikit-learn callback
infrastructure.

First a few imports and some random data for the rest of the script.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%

import numpy as np

from sklearn.base import BaseEstimator, clone
from sklearn.callback import CallbackSupportMixin, ProgressBar, with_callback_context
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.model_selection import check_cv
from sklearn.utils import check_random_state
from sklearn.utils.validation import check_is_fitted, validate_data

n_samples, n_features = 100, 4
rng = np.random.RandomState(42)
X = rng.rand(n_samples, n_features)


# %%
# Custom Estimator
# ----------------
# Here we demonstrate how to implement a custom estimator that supports callbacks. For
# the example, a simplified version of KMeans is presented. First, let's implement our
# `SimpleKmeans` estimator without the callback support.


class SimpleKmeans(BaseEstimator):
    def __init__(self, n_clusters=6, n_iter=100, random_state=None):
        self.n_clusters = n_clusters
        self.n_iter = n_iter
        self.random_state = random_state  # controls the centroids' initialization

    def _compute_labels(self, X):
        # Get the index of the closest centroid for each point in X.
        return np.argmin(euclidean_distances(X, self.cluster_centers_), axis=1)

    def fit(self, X, y=None):
        # `y` is not used but we need to declare it to adhere to scikit-learn's
        # estimators fit convention.

        # Input validation is a goood practice in estimators, for more information about
        # it you can refer to
        # https://scikit-learn.org/stable/developers/develop.html#input-validation.
        X = validate_data(self, X)
        random_state = check_random_state(self.random_state)
        # Randomnly initialize the centroids.
        self.cluster_centers_ = random_state.rand(self.n_clusters, X.shape[1])

        for i in range(self.n_iter):
            # The fit iterations consist in getting the cluster label of each data point
            # according to their closest centroid, and then updating the centroids as
            # the center of each cluster.
            labels = self._compute_labels(X)

            for k in range(self.n_clusters):
                # For each centroid, if its cluster is not empty, its coordinates are
                # updated with the coordinates of the cluster centers.
                if (labels == k).any():
                    self.cluster_centers_[k] = X[labels == k].mean(axis=0)

        return self

    def predict(self, X):
        check_is_fitted(self)
        return self._compute_labels(X)

    def transform(self, X):
        check_is_fitted(self)
        return euclidean_distances(X, self.cluster_centers_)


# %%
# Now let's add all the elements necessary to support callbacks.


# First things first, the estimator must inherit from the `CallbackSupportMixin` class.
class SimpleKMeans(CallbackSupportMixin, BaseEstimator):
    def __init__(self, n_clusters=6, n_iter=100, random_state=None):
        self.n_clusters = n_clusters
        self.n_iter = n_iter
        self.random_state = random_state

    def _compute_labels(self, X):
        return np.argmin(euclidean_distances(X, self.cluster_centers_), axis=1)

    # Then the `fit` function must be decorated with the `with_callback_context`
    # decorator, which will take care of the proper tear down of callbacks.
    @with_callback_context
    def fit(self, X, y=None):
        X = validate_data(self, X)
        random_state = check_random_state(self.random_state)
        # The `CallbackContext` object must be instantiated with the
        # `_init_callback_context` method provided by the mixin.
        callback_ctx = self._init_callback_context(
            task_name="fit", task_id=0, max_subtasks=self.n_iter
        )
        # Then the callback context's `eval_on_fit_begin` method must be called. It will
        # trigger all the callbacks' `on_fit_begin` methods.
        callback_ctx.eval_on_fit_begin(estimator=self)

        self.cluster_centers_ = random_state.rand(self.n_clusters, X.shape[1])

        for i in range(self.n_iter):
            # For each of the fit task (here each iteration of the loop), a subcontext
            # must be created with the callback context's `subcontext` method.
            subcontext = callback_ctx.subcontext(task_id=i, task_name="fit iteration")

            labels = self._compute_labels(X)

            for k in range(self.n_clusters):
                if (labels == k).any():
                    self.cluster_centers_[k] = X[labels == k].mean(axis=0)

            # After each task, the `eval_on_fit_task_end` method of its callback context
            # must be called. It will trigger all the callbacks' `on_fit_task_end`
            # methods. Extra `kwargs` can be passed to provide extra contextual tools
            # for the callbacks, such as the `reconstruction_attributes` function which
            # returns the necessary arguments to generate an estimator instance ready to
            # predict, as if the fit process just stopped at this step. See the
            # following note for more details on these `kwargs`.
            if subcontext.eval_on_fit_task_end(
                estimator=self,
                data={"X_train": X, "y_train": y},
                reconstruction_attributes=lambda: {
                    "centroids_": self.cluster_centers_,
                },
            ):
                # The `eval_on_fit_task_end` method returns a boolean, which will be set
                # to True if any of the callbacks' `on_fit_task_end` methods return
                # True. This allows the callbacks to interrupt the `fit` process, for
                # example to implement early stopping. Thus the `eval_on_fit_task_end`
                # method should be used in an `if` / `break` block to enable such
                # interruptions.
                break

        # The callback context's `eval_on_fit_end` method does not need to be called
        # here as it is called automatically in the decorator, after fit finishes, even
        # if it crashes. Thus the callbacks will call their `on_fit_end` methods even if
        # `fit` does not complete correctly.
        return self

    def predict(self, X):
        check_is_fitted(self)
        return self._compute_labels(X)

    def transform(self, X):
        check_is_fitted(self)
        return euclidean_distances(X, self.cluster_centers_)


# %%
# .. note::
#
#     See :func:`~sklearn.callback.CallbackContext.eval_on_fit_task_end` for a list of
#     the possible keys of the ``kwargs`` to provide to ``eval_on_fit_task_end``. These
#     ``kwargs`` are optional, but an estimator should provide all the ones it is
#     capable of producing to be compatible with a maximum number of callbacks.

# %%
# Registering callbacks to the custom estimator
# -------------------
# Now the ``SimpleKmeans`` estimator can be used with callbacks, for example with the
# :class:`~callback.ProgressBar` callback to monitor progress.

estimator = SimpleKMeans(random_state=rng)
callback = ProgressBar()
estimator.set_callbacks(callback)
estimator.fit(X)

# %%
# Custom meta-estimator
# ---------------------
# Now we will demonstrate how to implement a custom meta-estimator that supports
# callbacks. For the example, we will implement a simplified version of a grid search,
# where only a list of parameter combinations is provided and searched through instead
# of a grid. Again, let's start with the implementation without the callback support.


class SimpleGridSearch(BaseEstimator):
    def __init__(self, estimator, param_list, cv, score_func):
        # the estimator to evaluate
        self.estimator = estimator
        # the list of parameter combinations to iterate over
        self.param_list = param_list
        # the number of splits for the CV, or a CV splitter instance
        self.cv = cv
        # the scoring function
        self.score_func = score_func

    def fit(self, X, y=None):
        # We make a cross-validator instance to evaluate each parameter combination on
        # multiple folds.
        cv = check_cv(self.cv)
        # This `cv_results_` attribute will hold the score values for each parameter
        # combination and fold.
        self.cv_results_ = []

        # We iterate on the parameter combinations and the folds, computing a score
        # value for each param combination and fold.
        for i, params in enumerate(self.param_list):
            for j, (train_idx, test_idx) in enumerate(cv.split(X)):
                # An estimator is initialized with the parameter combination.
                estimator = clone(self.estimator).set_params(**params)
                # The split of the current fold is applied to the data.
                X_train, X_test = X[train_idx], X[test_idx]
                y_train, y_test = (
                    (y[train_idx], y[test_idx]) if y is not None else (None, None)
                )
                # The estimator is fitted.
                estimator.fit(X_train, y_train)
                # Its score is computed.
                score = self.score_func(estimator, X_test, y_test)
                # The results are aggregated as a tuple in the attribute.
                self.cv_results_.append((params, f"split_{j}", score))

        return self


# %%
# Now let's update the class to support callbacks.


# The class must again inherit from `CallbackSupportMixin`.
class SimpleGridSearch(CallbackSupportMixin, BaseEstimator):  # noqa: F811
    def __init__(self, estimator, param_list, cv, score_func):
        self.estimator = estimator
        self.param_list = param_list
        self.cv = cv
        self.score_func = score_func

    # The `fit` method must also be decorated.
    @with_callback_context
    def fit(self, X, y=None):
        cv = check_cv(self.cv)
        # The callback context must also be instantiated.
        callback_ctx = self._init_callback_context(
            task_name="fit", task_id=0, max_subtasks=len(self.param_list)
        )
        # The `eval_on_fit_begin` method must again be called.
        callback_ctx.eval_on_fit_begin(estimator=self)

        self.cv_results_ = []

        # The tasks of the `fit` function are here the iterations on two levels of
        # nested loops. Therefore, two levels of subcontexts must be used. Each level
        # will need its subcontext to call its `eval_on_fit_task_end` method. These
        # methods can be used at any level to enable interruptions through a `break`.
        # Here it does not make much sense to allow stopping between folds inside the
        # inner loop, so it is only implemented on the outer loop, but this is only a
        # design choice.

        for i, params in enumerate(self.param_list):
            # A subcontext for the first level of `fit` iterations must be created.
            outer_subcontext = callback_ctx.subcontext(
                task_name="param iteration", task_id=i, max_subtasks=cv.get_n_splits()
            )
            for j, (train_idx, test_idx) in enumerate(cv.split(X)):
                # This time a second level of `fit` iterations is also used, a second
                # level of subcontext must then be used.
                inner_subcontext = outer_subcontext.subcontext(
                    task_name=f"split {j}", task_id=j
                )

                estimator = clone(self.estimator).set_params(**params)
                # Since a sub-estimator is used, the callbacks must be propagated to
                # that estimator with the `propagate_callbacks` method. Note that only
                # the callbacks following the `AutoPropagatedCallback` protocol will be
                # propagated.
                inner_subcontext.propagate_callbacks(sub_estimator=estimator)

                X_train, X_test = X[train_idx], X[test_idx]
                y_train, y_test = (
                    (y[train_idx], y[test_idx]) if y is not None else (None, None)
                )
                estimator.fit(X_train, y_train)
                score = self.score_func(estimator, X_test, y_test)
                self.cv_results_.append((params, f"split_{j}", score))

                # The inner subcontext's `eval_on_fit_task_end` must be called.
                inner_subcontext.eval_on_fit_task_end(
                    estimator=self,
                    data={
                        "X_train": X_train,
                        "y_train": y_train,
                        "X_val": X_test,
                        "y_val": y_test,
                    },
                )

            # The outer subcontext's `eval_on_fit_task_end` must be called as well and
            # is used with an `if` / `break` to enable interruptions.
            if outer_subcontext.eval_on_fit_task_end(
                estimator=self,
                data={"X_train": X, "y_train": y, "X_val": None, "y_val": None},
            ):
                break

        # Again, the callback context's `eval_on_fit_end` is taken care of automatically
        # in the decorator.
        return self


# %%
# The main difference with a simple estimator is that the callbacks must be propagated
# to the sub-estimators through the corresponding callback subcontext's
# :meth:`~callback._callback_context.CallbackContext.propagate_callbacks` method.


# %%
# Registering callbacks to the meta-estimator
# -------------------------------------------
# Callbacks are registered to a meta-estimator the same way as to regular estimators.
# The callbacks which respect the :class:`~callback.AutoPropagatedCallback` protocol
# (such as :class:`~callback.ProgressBar`) will be propagated to the sub-estimators.


param_list = [{"n_clusters": 5, "n_iter": 20}, {"n_clusters": 4, "n_iter": 50}]


def score_func(estimator, X, y=None):
    return np.sum(estimator.transform(X).min(axis=1))


sub_estimator = SimpleKMeans(random_state=rng)
meta_estimator = SimpleGridSearch(
    estimator=sub_estimator, param_list=param_list, cv=4, score_func=score_func
)
callback = ProgressBar(max_estimator_depth=2)
meta_estimator.set_callbacks(callback)
meta_estimator.fit(X)
