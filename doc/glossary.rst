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
    common tests
        TODO
    cross validation splitter
        TODO
    double-underscore notation
        When specifying parameter names for nested estimators, ``__`` may be
        used to separate between parent and child.
        See :term:`parameter`.
    estimator
        TODO
    estimator instance
        We sometimes use this terminology to distinguish an estimator class
        from a constructed instance. For example, in the following, ``cls``
        is an estimator class, while ``est1`` and ``est2`` are instances::

            cls = RandomForestClassifier
            est1 = cls()
            est2 = RandomForestClassifier()
    feature
        TODO
    fitted
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
    predictor
        TODO
    sample
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

.. glossary::

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
        mention validation
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

.. glossary::

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

.. glossary::

    ``feature_importances_``
        TODO
    ``labels_``
        TODO

Sample properties
=================

.. glossary::

    ``groups``
        TODO
    ``sample_weight``
        TODO
