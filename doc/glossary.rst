.. _glossary:

=========================================
Glossary of Common Terms and API Elements
=========================================

This glossary hopes to definitively represent the tacit and explicit
conventions applied in scikit-learn and its API, while providing a reference
for users and contributors. It aims to describe the concepts and either detail
their corresponding API or link to other relevant parts of the documentation
which do so.

General Concepts
================

.. glossary::

    array-like
        TODO

    attribute
        TODO

    classifier
        TODO
        Mention :func:`~sklearn.base.is_classifier`.

    clone
        To copy an :term:`estimator instance` and create a new one with
        identical :term:`parameters <parameter>`, but without any fitted
        :term:`attributes <attribute>`, using :func:`~skelarn.base.clone`.

        When ``fit`` is called, a :term:`meta-estimator` usually clones
        a wrapped estimator instance before fitting the cloned instance.
        (Exceptions, for legacy reasons, include
        :class:`~sklearn.pipeline.Pipeline` and
        :class:`~sklearn.pipeline.FeatureUnion`.)

    common tests
        TODO

    convergence
        TODO mention :class:`sklearn.exceptions.ConvergenceWarning`

    deprecated
    deprecation
        TODO

    double underscore notation
        When specifying parameter names for nested estimators, ``__`` may be
        used to separate between parent and child.
        See :term:`parameter`.

    duck typing
        TODO
        Note that ``getattr`` should be preferred to ``hasattr`` since
        ``hasattr`` can be expensive, particularly for some model attributes.

    estimator instance
        We sometimes use this terminology to distinguish an :term:`estimator`
        class from a constructed instance. For example, in the following,
        ``cls`` is an estimator class, while ``est1`` and ``est2`` are
        instances::

            cls = RandomForestClassifier
            est1 = cls()
            est2 = RandomForestClassifier()

    evaluation metric
        TODO

    feature
        TODO

    fitted
        TODO

    missing data
        TODO

    parameter
        We mostly use *parameter* to refer to the aspects of an estimator that can be
        specified in its construction. For example, ``max_depth`` and ``random_state``
        are parameters of :class:`RandomForestClassifier`.

        We do not use _parameters_ in the statistical sense, where parameters
        are values that specify a model and can be estimated from data. In this
        sense, what we call parameters might be what statisticians call
        hyperparameters to the model: decisions about model structure that are
        often not directly learnt from data.  However, our parameters are also
        used to prescribe modeling operations that do not affect the learnt
        model, such as :term:`n_jobs` for controlling parallelism.

        When talking about the parameters of a :term:`meta-estimator`, we may
        also be including the parameters of the estimators wrapped by the
        meta-estimator.  Ordinarily, these nested parameters are denoted by
        using a double-underscore (``__``) to separate between the
        estimator-as-parameter and its parameter.  Thus
        ``BaggingClassifier(base_estimator=DecisionTreeClassifier(max_depth=3))``
        has a deep paramter ``base_estimator__max_depth`` with value ``3``.

        The list of parameters and their current values can be retrieved from
        an :term:`estimator instance` using its :term:`get_params` method.

        Between construction and fitting, parameters may be modified using
        :term:`set_params`.  To enable this, parameters are not ordinarily
        validated when the estimator is constructed, or when each parameter is
        set. Parameter validation is performed when :term:`fit` is called.

    pairwise metric
        TODO

    sample
        We usually use this terms as a noun to indicate a single instance or
        feature vector.  Thus ``n_samples`` indicates the number of instances
        in a dataset.

    sample property
        TODO

    scikit-learn-contrib
        TODO

    scorer
        TODO
        See also :term:`evaluation metric`.

    target
        TODO

    unlabeled data
        TODO

Class APIs and Estimator Types
==============================

.. glossary::

    cross validation splitter
        TODO

    estimator
        TODO

    meta-estimator
        TODO
        Mention duck typing.

    outlier detector
        TODO

    predictor
        TODO

    regressor
        TODO
        Mention :func:`~sklearn.base.is_regressor`.

    transformer
        TODO

    vectorizer
        TODO

Target Types
============

.. glossary::

    binary
        TODO

    continuous
        TODO

    multiclass
        TODO

    multilabel
        TODO

    multi-output continuous
        TODO

    multi-output multiclass
        TODO

Methods
=======

.. glossary::

    ``decision_function``

        classifier
            TODO

        outlier detector
            TODO

    ``get_feature_names``
        TODO

    ``get_n_splits``
        TODO

    ``get_params``
        TODO

    ``fit_predict``
        TODO

    ``fit_transform``
        TODO

    ``fit``
        TODO
        mention validation

    ``partial_fit``
        TODO

    ``predict_log_proba``
        TODO

    ``predict_proba``
        TODO

    ``predict``
        TODO
        Mention ``return_std``

    ``score``
        TODO

    ``score_samples``
        TODO

    ``set_params``
        TODO

    ``split``
        TODO

Parameters
==========

See concept :term:`parameter`.

.. glossary::

    ``cv``
        TODO

    ``max_iter``
        TODO

    ``n_iter``
        TODO

    ``n_jobs``
        TODO

    ``random_state``
        TODO

    ``scoring``
        TODO

    ``verbose``
        TODO

    ``warm_start``
        TODO
        See also :term:`partial_fit`.

Attributes
==========

See concept :term:`attribute`.

.. glossary::

    ``n_iter_``
        TODO

    ``feature_importances_``
        TODO

    ``labels_``
        TODO

Sample properties
=================

See concept :term:`sample property`.

.. glossary::

    ``groups``
        TODO

    ``sample_weight``
        TODO
