.. _glossary:

=========================================
Glossary of Common Terms and API Elements
=========================================

This glossary hopes to definitively represent the tacit and explicit standards
applied in scikit-learn and its API, while providing a reference for users and
contributors. It aims to describe the concepts and either detail their
corresponding API or link to other relevant parts of the documentation which do
so.

Concepts
========

.. glossary::

    attribute (for representing a learnt model)
        TODO
    clone
        TODO
    cross validation splitter
        TODO
    estimator
        TODO
    meta-estimator
        TODO
    metric
        TODO
    multiclass
        TODO
    multilabel
        TODO
    multi-output
        TODO
    parameter
        We mostly use parameter to refer to the aspects of an estimator that can be
        specified in construction. For example, ``min_depth`` and ``random_state``
        are parameter of ``RandomForestClassifier``.

        We do not use _parameters_ in the statistical sense, where parameters
        are values that specify a model and can be estimated from data. In this
        sense, what we call parameters might be what statisticians call
        hyperparameters to the model: decisions about model structure that are
        often not directly learnt from data.
    pairwise metric
        TODO
    predictor
        TODO
    sample property
        TODO
    scikit-learn-contrib
        TODO
    scorer
        TODO
        See also :term:`metric`.
    target
        TODO
    transformer
        TODO
    vectorizer
        TODO

Methods
=======

    ``decision_function``
        TODO
    ``get_params``
        TODO
    ``fit_predict``
        TODO
    ``fit_transform``
        TODO
    ``fit``
        TODO
    ``partial_fit``
        TODO
    ``predict_log_proba``
        TODO
    ``predict_proba``
        TODO
    ``predict``
        TODO
    ``set_params``
        TODO

Parameters
==========

    ``n_jobs``
        TODO
    ``random_state``
        TODO
    ``scoring``
        TODO
    ``warm_start``
        TODO
        See also :term:`partial_fit`.

Attributes
==========

    ``feature_importances_``
        TODO
    ``labels_``
        TODO

Sample properties
=================

    ``groups``
        TODO
    ``sample_weight``
        TODO
