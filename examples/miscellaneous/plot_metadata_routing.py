"""
================
Metadata Routing
================

.. currentmodule:: sklearn

This document shows how you can use the :ref:`metadata routing mechanism
<metadata_routing>` in scikit-learn to route metadata to the estimators,
scorers, and CV splitters consuming them.

To better understand the following document, we need to introduce two concepts:
routers and consumers. A router is an object which forwards some given data and
metadata to other objects. In most cases, a router is a :term:`meta-estimator`,
i.e. an estimator which takes another estimator as a parameter. A function such
as :func:`sklearn.model_selection.cross_validate` which takes an estimator as a
parameter and forwards data and metadata, is also a router.

A consumer, on the other hand, is an object which accepts and uses some given
metadata. For instance, an estimator taking into account ``sample_weight`` in
its :term:`fit` method is a consumer of ``sample_weight``.

It is possible for an object to be both a router and a consumer. For instance,
a meta-estimator may take into account ``sample_weight`` in certain
calculations, but it may also route it to the underlying estimator.

First a few imports and some random data for the rest of the script.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%

import warnings
from pprint import pprint

import numpy as np

from sklearn import set_config
from sklearn.base import (
    BaseEstimator,
    ClassifierMixin,
    MetaEstimatorMixin,
    RegressorMixin,
    TransformerMixin,
    clone,
)
from sklearn.linear_model import LinearRegression
from sklearn.utils import metadata_routing
from sklearn.utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    get_routing_for_object,
    process_routing,
)
from sklearn.utils.validation import check_is_fitted

n_samples, n_features = 100, 4
rng = np.random.RandomState(42)
X = rng.rand(n_samples, n_features)
y = rng.randint(0, 2, size=n_samples)
my_groups = rng.randint(0, 10, size=n_samples)
my_weights = rng.rand(n_samples)
my_other_weights = rng.rand(n_samples)

# %%
# Metadata routing is only available if explicitly enabled:
set_config(enable_metadata_routing=True)


# %%
# This utility function is a dummy to check if a metadata is passed:
def check_metadata(obj, **kwargs):
    for key, value in kwargs.items():
        if value is not None:
            print(
                f"Received {key} of length = {len(value)} in {obj.__class__.__name__}."
            )
        else:
            print(f"{key} is None in {obj.__class__.__name__}.")


# %%
# A utility function to nicely print the routing information of an object:
def print_routing(obj):
    pprint(obj.get_metadata_routing()._serialize())


# %%
# Consuming Estimator
# -------------------
# Here we demonstrate how an estimator can expose the required API to support
# metadata routing as a consumer. Imagine a simple classifier accepting
# ``sample_weight`` as a metadata on its ``fit`` and ``groups`` in its
# ``predict`` method:


class ExampleClassifier(ClassifierMixin, BaseEstimator):
    def fit(self, X, y, sample_weight=None):
        check_metadata(self, sample_weight=sample_weight)
        # all classifiers need to expose a classes_ attribute once they're fit.
        self.classes_ = np.array([0, 1])
        return self

    def predict(self, X, groups=None):
        check_metadata(self, groups=groups)
        # return a constant value of 1, not a very smart classifier!
        return np.ones(len(X))


# %%
# The above estimator now has all it needs to consume metadata. This is
# accomplished by some magic done in :class:`~base.BaseEstimator`. There are
# now three methods exposed by the above class: ``set_fit_request``,
# ``set_predict_request``, and ``get_metadata_routing``. There is also a
# ``set_score_request`` for ``sample_weight`` which is present since
# :class:`~base.ClassifierMixin` implements a ``score`` method accepting
# ``sample_weight``. The same applies to regressors which inherit from
# :class:`~base.RegressorMixin`.
#
# By default, no metadata is requested, which we can see as:

print_routing(ExampleClassifier())

# %%
# The above output means that ``sample_weight`` and ``groups`` are not
# requested by `ExampleClassifier`, and if a router is given those metadata, it
# should raise an error, since the user has not explicitly set whether they are
# required or not. The same is true for ``sample_weight`` in the ``score``
# method, which is inherited from :class:`~base.ClassifierMixin`. In order to
# explicitly set request values for those metadata, we can use these methods:

est = (
    ExampleClassifier()
    .set_fit_request(sample_weight=False)
    .set_predict_request(groups=True)
    .set_score_request(sample_weight=False)
)
print_routing(est)

# %%
# .. note ::
#     Please note that as long as the above estimator is not used in a
#     meta-estimator, the user does not need to set any requests for the
#     metadata and the set values are ignored, since a consumer does not
#     validate or route given metadata. A simple usage of the above estimator
#     would work as expected.

est = ExampleClassifier()
est.fit(X, y, sample_weight=my_weights)
est.predict(X[:3, :], groups=my_groups)

# %%
# Routing Meta-Estimator
# ----------------------
# Now, we show how to design a meta-estimator to be a router. As a simplified
# example, here is a meta-estimator, which doesn't do much other than routing
# the metadata.


class MetaClassifier(MetaEstimatorMixin, ClassifierMixin, BaseEstimator):
    def __init__(self, estimator):
        self.estimator = estimator

    def get_metadata_routing(self):
        # This method defines the routing for this meta-estimator.
        # In order to do so, a `MetadataRouter` instance is created, and the
        # routing is added to it. More explanations follow below.
        router = MetadataRouter(owner=self.__class__.__name__).add(
            estimator=self.estimator,
            method_mapping=MethodMapping()
            .add(caller="fit", callee="fit")
            .add(caller="predict", callee="predict")
            .add(caller="score", callee="score"),
        )
        return router

    def fit(self, X, y, **fit_params):
        # `get_routing_for_object` returns a copy of the `MetadataRouter`
        # constructed by the above `get_metadata_routing` method, that is
        # internally called.
        request_router = get_routing_for_object(self)
        # Meta-estimators are responsible for validating the given metadata.
        # `method` refers to the parent's method, i.e. `fit` in this example.
        request_router.validate_metadata(params=fit_params, method="fit")
        # `MetadataRouter.route_params` maps the given metadata to the metadata
        # required by the underlying estimator based on the routing information
        # defined by the MetadataRouter. The output of type `Bunch` has a key
        # for each consuming object and those hold keys for their consuming
        # methods, which then contain key for the metadata which should be
        # routed to them.
        routed_params = request_router.route_params(params=fit_params, caller="fit")

        # A sub-estimator is fitted and its classes are attributed to the
        # meta-estimator.
        self.estimator_ = clone(self.estimator).fit(X, y, **routed_params.estimator.fit)
        self.classes_ = self.estimator_.classes_
        return self

    def predict(self, X, **predict_params):
        check_is_fitted(self)
        # As in `fit`, we get a copy of the object's MetadataRouter,
        request_router = get_routing_for_object(self)
        # then we validate the given metadata,
        request_router.validate_metadata(params=predict_params, method="predict")
        # and then prepare the input to the underlying `predict` method.
        routed_params = request_router.route_params(
            params=predict_params, caller="predict"
        )
        return self.estimator_.predict(X, **routed_params.estimator.predict)


# %%
# Let's break down different parts of the above code.
#
# First, the :meth:`~utils.metadata_routing.get_routing_for_object` takes our
# meta-estimator (``self``) and returns a
# :class:`~utils.metadata_routing.MetadataRouter` or, a
# :class:`~utils.metadata_routing.MetadataRequest` if the object is a consumer,
# based on the output of the estimator's ``get_metadata_routing`` method.
#
# Then in each method, we use the ``route_params`` method to construct a
# dictionary of the form ``{"object_name": {"method_name": {"metadata":
# value}}}`` to pass to the underlying estimator's method. The ``object_name``
# (``estimator`` in the above ``routed_params.estimator.fit`` example) is the
# same as the one added in the ``get_metadata_routing``. ``validate_metadata``
# makes sure all given metadata are requested to avoid silent bugs.
#
# Next, we illustrate the different behaviors and notably the type of errors
# raised.

meta_est = MetaClassifier(
    estimator=ExampleClassifier().set_fit_request(sample_weight=True)
)
meta_est.fit(X, y, sample_weight=my_weights)

# %%
# Note that the above example is calling our utility function
# `check_metadata()` via the `ExampleClassifier`. It checks that
# ``sample_weight`` is correctly passed to it. If it is not, like in the
# following example, it would print that ``sample_weight`` is ``None``:

meta_est.fit(X, y)

# %%
# If we pass an unknown metadata, an error is raised:
try:
    meta_est.fit(X, y, test=my_weights)
except TypeError as e:
    print(e)

# %%
# And if we pass a metadata which is not explicitly requested:
try:
    meta_est.fit(X, y, sample_weight=my_weights).predict(X, groups=my_groups)
except ValueError as e:
    print(e)

# %%
# Also, if we explicitly set it as not requested, but it is provided:
meta_est = MetaClassifier(
    estimator=ExampleClassifier()
    .set_fit_request(sample_weight=True)
    .set_predict_request(groups=False)
)
try:
    meta_est.fit(X, y, sample_weight=my_weights).predict(X[:3, :], groups=my_groups)
except TypeError as e:
    print(e)

# %%
# Another concept to introduce is **aliased metadata**. This is when an
# estimator requests a metadata with a different variable name than the default
# variable name. For instance, in a setting where there are two estimators in a
# pipeline, one could request ``sample_weight1`` and the other
# ``sample_weight2``. Note that this doesn't change what the estimator expects,
# it only tells the meta-estimator how to map the provided metadata to what is
# required. Here's an example, where we pass ``aliased_sample_weight`` to the
# meta-estimator, but the meta-estimator understands that
# ``aliased_sample_weight`` is an alias for ``sample_weight``, and passes it as
# ``sample_weight`` to the underlying estimator:
meta_est = MetaClassifier(
    estimator=ExampleClassifier().set_fit_request(sample_weight="aliased_sample_weight")
)
meta_est.fit(X, y, aliased_sample_weight=my_weights)

# %%
# Passing ``sample_weight`` here will fail since it is requested with an
# alias and ``sample_weight`` with that name is not requested:
try:
    meta_est.fit(X, y, sample_weight=my_weights)
except TypeError as e:
    print(e)

# %%
# This leads us to the ``get_metadata_routing``. The way routing works in
# scikit-learn is that consumers request what they need, and routers pass that
# along. Additionally, a router exposes what it requires itself so that it can
# be used inside another router, e.g. a pipeline inside a grid search object.
# The output of the ``get_metadata_routing`` which is a dictionary
# representation of a :class:`~utils.metadata_routing.MetadataRouter`, includes
# the complete tree of requested metadata by all nested objects and their
# corresponding method routings, i.e. which method of a sub-estimator is used
# in which method of a meta-estimator:

print_routing(meta_est)

# %%
# As you can see, the only metadata requested for method ``fit`` is
# ``"sample_weight"`` with ``"aliased_sample_weight"`` as the alias. The
# ``~utils.metadata_routing.MetadataRouter`` class enables us to easily create
# the routing object which would create the output we need for our
# ``get_metadata_routing``.
#
# In order to understand how aliases work in meta-estimators, imagine our
# meta-estimator inside another one:

meta_meta_est = MetaClassifier(estimator=meta_est).fit(
    X, y, aliased_sample_weight=my_weights
)

# %%
# In the above example, this is how the ``fit`` method of `meta_meta_est`
# will call their sub-estimator's ``fit`` methods::
#
#     # user feeds `my_weights` as `aliased_sample_weight` into `meta_meta_est`:
#     meta_meta_est.fit(X, y, aliased_sample_weight=my_weights):
#         ...
#
#         # the first sub-estimator (`meta_est`) expects `aliased_sample_weight`
#         self.estimator_.fit(X, y, aliased_sample_weight=aliased_sample_weight):
#             ...
#
#             # the second sub-estimator (`est`) expects `sample_weight`
#             self.estimator_.fit(X, y, sample_weight=aliased_sample_weight):
#                 ...

# %%
# Consuming and routing Meta-Estimator
# ------------------------------------
# For a slightly more complex example, consider a meta-estimator that routes
# metadata to an underlying estimator as before, but it also uses some metadata
# in its own methods. This meta-estimator is a consumer and a router at the
# same time. Implementing one is very similar to what we had before, but with a
# few tweaks.


class RouterConsumerClassifier(MetaEstimatorMixin, ClassifierMixin, BaseEstimator):
    def __init__(self, estimator):
        self.estimator = estimator

    def get_metadata_routing(self):
        router = (
            MetadataRouter(owner=self.__class__.__name__)
            # defining metadata routing request values for usage in the meta-estimator
            .add_self_request(self)
            # defining metadata routing request values for usage in the sub-estimator
            .add(
                estimator=self.estimator,
                method_mapping=MethodMapping()
                .add(caller="fit", callee="fit")
                .add(caller="predict", callee="predict")
                .add(caller="score", callee="score"),
            )
        )
        return router

    # Since `sample_weight` is used and consumed here, it should be defined as
    # an explicit argument in the method's signature. All other metadata which
    # are only routed, will be passed as `**fit_params`:
    def fit(self, X, y, sample_weight, **fit_params):
        if self.estimator is None:
            raise ValueError("estimator cannot be None!")

        check_metadata(self, sample_weight=sample_weight)

        # We add `sample_weight` to the `fit_params` dictionary.
        if sample_weight is not None:
            fit_params["sample_weight"] = sample_weight

        request_router = get_routing_for_object(self)
        request_router.validate_metadata(params=fit_params, method="fit")
        routed_params = request_router.route_params(params=fit_params, caller="fit")
        self.estimator_ = clone(self.estimator).fit(X, y, **routed_params.estimator.fit)
        self.classes_ = self.estimator_.classes_
        return self

    def predict(self, X, **predict_params):
        check_is_fitted(self)
        # As in `fit`, we get a copy of the object's MetadataRouter,
        request_router = get_routing_for_object(self)
        # we validate the given metadata,
        request_router.validate_metadata(params=predict_params, method="predict")
        # and then prepare the input to the underlying ``predict`` method.
        routed_params = request_router.route_params(
            params=predict_params, caller="predict"
        )
        return self.estimator_.predict(X, **routed_params.estimator.predict)


# %%
# The key parts where the above meta-estimator differs from our previous
# meta-estimator is accepting ``sample_weight`` explicitly in ``fit`` and
# including it in ``fit_params``. Since ``sample_weight`` is an explicit
# argument, we can be sure that ``set_fit_request(sample_weight=...)`` is
# present for this method. The meta-estimator is both a consumer, as well as a
# router of ``sample_weight``.
#
# In ``get_metadata_routing``, we add ``self`` to the routing using
# ``add_self_request`` to indicate this estimator is consuming
# ``sample_weight`` as well as being a router; which also adds a
# ``$self_request`` key to the routing info as illustrated below. Now let's
# look at some examples:

# %%
# - No metadata requested
meta_est = RouterConsumerClassifier(estimator=ExampleClassifier())
print_routing(meta_est)


# %%
# - ``sample_weight`` requested by sub-estimator
meta_est = RouterConsumerClassifier(
    estimator=ExampleClassifier().set_fit_request(sample_weight=True)
)
print_routing(meta_est)

# %%
# - ``sample_weight`` requested by meta-estimator
meta_est = RouterConsumerClassifier(estimator=ExampleClassifier()).set_fit_request(
    sample_weight=True
)
print_routing(meta_est)

# %%
# Note the difference in the requested metadata representations above.
#
# - We can also alias the metadata to pass different values to the fit methods
#   of the meta- and the sub-estimator:

meta_est = RouterConsumerClassifier(
    estimator=ExampleClassifier().set_fit_request(sample_weight="clf_sample_weight"),
).set_fit_request(sample_weight="meta_clf_sample_weight")
print_routing(meta_est)

# %%
# However, ``fit`` of the meta-estimator only needs the alias for the
# sub-estimator and addresses their own sample weight as `sample_weight`, since
# it doesn't validate and route its own required metadata:
meta_est.fit(X, y, sample_weight=my_weights, clf_sample_weight=my_other_weights)

# %%
# - Alias only on the sub-estimator:
#
# This is useful when we don't want the meta-estimator to use the metadata, but
# the sub-estimator should.
meta_est = RouterConsumerClassifier(
    estimator=ExampleClassifier().set_fit_request(sample_weight="aliased_sample_weight")
)
print_routing(meta_est)
# %%
# The meta-estimator cannot use `aliased_sample_weight`, because it expects
# it passed as `sample_weight`. This would apply even if
# `set_fit_request(sample_weight=True)` was set on it.

# %%
# Simple Pipeline
# ---------------
# A slightly more complicated use-case is a meta-estimator resembling a
# :class:`~pipeline.Pipeline`. Here is a meta-estimator, which accepts a
# transformer and a classifier. When calling its `fit` method, it applies the
# transformer's `fit` and `transform` before running the classifier on the
# transformed data. Upon `predict`, it applies the transformer's `transform`
# before predicting with the classifier's `predict` method on the transformed
# new data.


class SimplePipeline(ClassifierMixin, BaseEstimator):
    def __init__(self, transformer, classifier):
        self.transformer = transformer
        self.classifier = classifier

    def get_metadata_routing(self):
        router = (
            MetadataRouter(owner=self.__class__.__name__)
            # We add the routing for the transformer.
            .add(
                transformer=self.transformer,
                method_mapping=MethodMapping()
                # The metadata is routed such that it retraces how
                # `SimplePipeline` internally calls the transformer's `fit` and
                # `transform` methods in its own methods (`fit` and `predict`).
                .add(caller="fit", callee="fit")
                .add(caller="fit", callee="transform")
                .add(caller="predict", callee="transform"),
            )
            # We add the routing for the classifier.
            .add(
                classifier=self.classifier,
                method_mapping=MethodMapping()
                .add(caller="fit", callee="fit")
                .add(caller="predict", callee="predict"),
            )
        )
        return router

    def fit(self, X, y, **fit_params):
        routed_params = process_routing(self, "fit", **fit_params)

        self.transformer_ = clone(self.transformer).fit(
            X, y, **routed_params.transformer.fit
        )
        X_transformed = self.transformer_.transform(
            X, **routed_params.transformer.transform
        )

        self.classifier_ = clone(self.classifier).fit(
            X_transformed, y, **routed_params.classifier.fit
        )
        return self

    def predict(self, X, **predict_params):
        routed_params = process_routing(self, "predict", **predict_params)

        X_transformed = self.transformer_.transform(
            X, **routed_params.transformer.transform
        )
        return self.classifier_.predict(
            X_transformed, **routed_params.classifier.predict
        )


# %%
# Note the usage of :class:`~utils.metadata_routing.MethodMapping` to
# declare which methods of the child estimator (callee) are used in which
# methods of the meta estimator (caller). As you can see, `SimplePipeline` uses
# the transformer's ``transform`` and ``fit`` methods in ``fit``, and its
# ``transform`` method in ``predict``, and that's what you see implemented in
# the routing structure of the pipeline class.
#
# Another difference in the above example with the previous ones is the usage
# of :func:`~utils.metadata_routing.process_routing`, which processes the input
# parameters, does the required validation, and returns the `routed_params`
# which we had created in previous examples. This reduces the boilerplate code
# a developer needs to write in each meta-estimator's method. Developers are
# strongly recommended to use this function unless there is a good reason
# against it.
#
# In order to test the above pipeline, let's add an example transformer.


class ExampleTransformer(TransformerMixin, BaseEstimator):
    def fit(self, X, y, sample_weight=None):
        check_metadata(self, sample_weight=sample_weight)
        return self

    def transform(self, X, groups=None):
        check_metadata(self, groups=groups)
        return X

    def fit_transform(self, X, y, sample_weight=None, groups=None):
        return self.fit(X, y, sample_weight).transform(X, groups)


# %%
# Note that in the above example, we have implemented ``fit_transform`` which
# calls ``fit`` and ``transform`` with the appropriate metadata. This is only
# required if ``transform`` accepts metadata, since the default ``fit_transform``
# implementation in :class:`~base.TransformerMixin` doesn't pass metadata to
# ``transform``.
#
# Now we can test our pipeline, and see if metadata is correctly passed around.
# This example uses our `SimplePipeline`, our `ExampleTransformer`, and our
# `RouterConsumerClassifier` which uses our `ExampleClassifier`.

pipe = SimplePipeline(
    transformer=ExampleTransformer()
    # we set transformer's fit to receive sample_weight
    .set_fit_request(sample_weight=True)
    # we set transformer's transform to receive groups
    .set_transform_request(groups=True),
    classifier=RouterConsumerClassifier(
        estimator=ExampleClassifier()
        # we want this sub-estimator to receive sample_weight in fit
        .set_fit_request(sample_weight=True)
        # but not groups in predict
        .set_predict_request(groups=False),
    )
    # and we want the meta-estimator to receive sample_weight as well
    .set_fit_request(sample_weight=True),
)
pipe.fit(X, y, sample_weight=my_weights, groups=my_groups).predict(
    X[:3], groups=my_groups
)

# %%
# Deprecation / Default Value Change
# ----------------------------------
# In this section we show how one should handle the case where a router becomes
# also a consumer, especially when it consumes the same metadata as its
# sub-estimator, or a consumer starts consuming a metadata which it wasn't in
# an older release. In this case, a warning should be raised for a while, to
# let users know the behavior is changed from previous versions.


class MetaRegressor(MetaEstimatorMixin, RegressorMixin, BaseEstimator):
    def __init__(self, estimator):
        self.estimator = estimator

    def fit(self, X, y, **fit_params):
        routed_params = process_routing(self, "fit", **fit_params)
        self.estimator_ = clone(self.estimator).fit(X, y, **routed_params.estimator.fit)

    def get_metadata_routing(self):
        router = MetadataRouter(owner=self.__class__.__name__).add(
            estimator=self.estimator,
            method_mapping=MethodMapping().add(caller="fit", callee="fit"),
        )
        return router


# %%
# As explained above, this is a valid usage if `my_weights` aren't supposed
# to be passed as `sample_weight` to `MetaRegressor`:

reg = MetaRegressor(estimator=LinearRegression().set_fit_request(sample_weight=True))
reg.fit(X, y, sample_weight=my_weights)


# %%
# Now imagine we further develop ``MetaRegressor`` and it now also *consumes*
# ``sample_weight``:


class WeightedMetaRegressor(MetaEstimatorMixin, RegressorMixin, BaseEstimator):
    # show warning to remind user to explicitly set the value with
    # `.set_{method}_request(sample_weight={boolean})`
    __metadata_request__fit = {"sample_weight": metadata_routing.WARN}

    def __init__(self, estimator):
        self.estimator = estimator

    def fit(self, X, y, sample_weight=None, **fit_params):
        routed_params = process_routing(
            self, "fit", sample_weight=sample_weight, **fit_params
        )
        check_metadata(self, sample_weight=sample_weight)
        self.estimator_ = clone(self.estimator).fit(X, y, **routed_params.estimator.fit)

    def get_metadata_routing(self):
        router = (
            MetadataRouter(owner=self.__class__.__name__)
            .add_self_request(self)
            .add(
                estimator=self.estimator,
                method_mapping=MethodMapping().add(caller="fit", callee="fit"),
            )
        )
        return router


# %%
# The above implementation is almost the same as ``MetaRegressor``, and
# because of the default request value defined in ``__metadata_request__fit``
# there is a warning raised when fitted.

with warnings.catch_warnings(record=True) as record:
    WeightedMetaRegressor(
        estimator=LinearRegression().set_fit_request(sample_weight=False)
    ).fit(X, y, sample_weight=my_weights)
for w in record:
    print(w.message)


# %%
# When an estimator consumes a metadata which it didn't consume before, the
# following pattern can be used to warn the users about it.


class ExampleRegressor(RegressorMixin, BaseEstimator):
    __metadata_request__fit = {"sample_weight": metadata_routing.WARN}

    def fit(self, X, y, sample_weight=None):
        check_metadata(self, sample_weight=sample_weight)
        return self

    def predict(self, X):
        return np.zeros(shape=(len(X)))


with warnings.catch_warnings(record=True) as record:
    MetaRegressor(estimator=ExampleRegressor()).fit(X, y, sample_weight=my_weights)
for w in record:
    print(w.message)

# %%
# At the end we disable the configuration flag for metadata routing:

set_config(enable_metadata_routing=False)

# %%
# Third Party Development and scikit-learn Dependency
# ---------------------------------------------------
#
# As seen above, information is communicated between classes using
# :class:`~utils.metadata_routing.MetadataRequest` and
# :class:`~utils.metadata_routing.MetadataRouter`. It is strongly not advised,
# but possible to vendor the tools related to metadata-routing if you strictly
# want to have a scikit-learn compatible estimator, without depending on the
# scikit-learn package. If all of the following conditions are met, you do NOT
# need to modify your code at all:
#
# - your estimator inherits from :class:`~base.BaseEstimator`
# - the parameters consumed by your estimator's methods, e.g. ``fit``, are
#   explicitly defined in the method's signature, as opposed to being
#   ``*args`` or ``*kwargs``.
# - your estimator does not route any metadata to the underlying objects, i.e.
#   it's not a *router*.
