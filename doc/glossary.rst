.. currentmodule:: sklearn

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
    attributes
        TODO

        Sufficient statistics for prediction/transformation. Diagnostics.

    classifier
        TODO
        Mention :func:`~base.is_classifier`.

    clone
        To copy an :term:`estimator instance` and create a new one with
        identical :term:`parameters`, but without any fitted
        :term:`attributes`, using :func:`~skelarn.base.clone`.

        When ``fit`` is called, a :term:`meta-estimator` usually clones
        a wrapped estimator instance before fitting the cloned instance.
        (Exceptions, for legacy reasons, include
        :class:`~pipeline.Pipeline` and
        :class:`~pipeline.FeatureUnion`.)

    common tests
        TODO

    convergence
        TODO mention :class:`exceptions.ConvergenceWarning`

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

    fitting
        Calling :term:`fit` on an estimator.

    fitted
        TODO

    missing data
        TODO

    parameter
    parameters
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
        validated or altered when the estimator is constructed, or when each
        parameter is set. Parameter validation is performed when :term:`fit` is
        called.

    pairwise metric
        TODO

    sample
    samples
        We usually use this terms as a noun to indicate a single instance or
        feature vector.  Thus ``n_samples`` indicates the number of instances
        in a dataset.

    sample property
    sample properties
        TODO

    scikit-learn-contrib
        TODO

    scorer
        TODO
        See :ref:`scoring_parameter`.
        See also :term:`evaluation metric`.

    sparse matrix
        TODO

    target
    targets
        TODO

    unlabeled data
        TODO

Class APIs and Estimator Types
==============================

.. glossary::

    clusterer
        TODO

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
        An :term:`estimator` which provides :term:`predict`.
        This encompasses :term:`classifier`, :term:`regressor`,
        :term:`outlier detector` and sometimes :term:`clusterer` (at least when
        they are inductive).  In scikit-learn, if an estimator is not a
        predictor, it is usually a :term:`transformer`.

    regressor
        TODO
        Mention :func:`~base.is_regressor`.

    transformer
        TODO

    vectorizer
        A :term:`tranformer` which takes input where each sample is not
        represented as an :term:`array-like` object of fixed length, and
        produces an `array-like` object for each sample (and thus a
        2-dimensional array-like for a set of samples).

There are further APIs specific to a small family of estimators, such as:

* :class:`neighbors.DistanceMetric`
* :class:`gaussian_process.kernels.Kernel`
* ``tree.Criterion``

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
        The ``fit`` method is provided on every estimator. It usually takes some
        :term`samples` ``X``, :term:`targets` ``y`` if the model is supervised,
        and potentially other :term:`sample properties` such as
        :term:`sample_weight`.  It should:

        * clear any prior :term:`attributes` stored on the estimator, unless
          :term:`warm_start` is used;
        * validate and interpret any :term:`parameters`, ideally raising an
          error if invalid;
        * validate the input data;
        * estimate and store model attributes from the estimated parameters and
          provided data; and
        * return the now :term:`fitted` estimator to facilitate method
          chaining.

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

These common parameter names, specifically used in estimator construction
(see concept :term:`parameter`) sometimes also appear as function and
non-estimator parameters with similar semantics.

.. glossary::

    ``class_weight``
        TODO

    ``cv``
        TODO

    ``max_iter``
        TODO

        Mention ConvergenceWarning

    ``memory``
        Some estimators make use of :class:`joblib.Memory` to
        store partial solutions during fitting. Thus when ``fit`` is called
        again, those partial solutions have been memoized and can be reused.

        Memory can be specified as a string with a path to a directory, or
        a :class:`joblib.Memory` instance can be used. In the latter case,
        ``Memory`` should be imported from ``sklearn.externals.joblib``, which
        may differ from the installed ``joblib`` package.

    ``n_jobs``
        This is used to specify how many concurrent processes/threads should be
        used for parallelized routines.  Scikit-learn uses one processor for
        its processing by default, although it also makes use of NumPy, which
        may be configured to use a threaded numerical processor library (like
        MKL).

        ``n_jobs`` is an int, specifying the maximum number of concurrently
        running jobs.  If set to -1, all CPUs are used. If 1 is given, no
        parallel computing code is used at all.  For n_jobs below -1, (n_cpus +
        1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are used.

        The use of ``n_jobs``-based parallelism in estimators varies:

        * Most often parallelism happens in :term:`fitting <fit>`, but
          sometimes parallelism happens in prediction (e.g. in random forests).
        * Some parallelism uses a multi-threading backend by default, some
          a multi-processing backend.  It is possible to override the default
          backend by using :func:`sklearn.externals.joblib.parallel.parallel_backend`.
        * Whether parallel processing is helpful at improving runtime depends
          on many factors, and it's usually a good idea to experiment rather
          than assuming that increasing the number of jobs is always a good
          thing.

        Nested uses of ``n_jobs``-based parallelism with the same backend will
        result in an exception.
        So ``GridSearchCV(OneVsRestClassifier(SVC(), n_jobs=2), n_jobs=2)``
        won't work.

        When ``n_jobs`` is not 1, the estimator being parallelized must be
        picklable.  This means, for instance, that lambdas cannot be used
        as estimator parameters.

    ``random_state``
        Whenever randomization is part of a Scikit-learn algorithm, a
        ``random_state`` parameter may be provided to control the random number
        generator used.  Note that the mere presence of ``random_state`` doesn't
        mean that randomization is always used, as it may be dependent on
        another parameter, e.g. ``shuffle``, being set.

        ``random_state``'s value may be:

            None (default)
                Use the global random state from :mod:`numpy.random`.

            An integer
                Use a new random number generator seeded by the given integer.
                To make a randomized algorithm deterministic (i.e. running it
                multiple times will produce the same result), an arbitrary
                integer ``random_state`` can be used. However, it may be
                worthwhile checking that your results are stable across a
                number of different distinct random seeds.

            A :class:`numpy.random.RandomState` instance
                Use the provided random state, only affecting other users
                of the same random state instance.

        :func:`utils.check_random_state` is used internally to validate the
        input ``random_state`` and return a :class:`~numpy.random.RandomState`
        instance.

    ``scoring``
        Specifies the score function to be maximized (usually by :ref:`cross
        validation <cross_validation>`), or -- in some cases -- multiple score
        functions to be reported. The score function can be a string accepted
        by :func:`metrics.get_scorer` or a callable :term:`scorer`, not to be
        confused with an :term:`evaluation metric`, as the latter have a more
        diverse API.  ``scoring`` may also be set to None, in which case the
        estimator's ``score`` method is used.  See :ref:`scoring_parameter` in
        the user guide.

        Where multiple metrics can be evaluated, ``scoring`` may be given
        either as a list of unique strings or a dict with names as keys and
        callables as values. Note that this does *not* specify which score
        function is to be maximised, and another parameter such as ``refit``
        may be used for this purpose.

    ``verbose``
        TODO

    ``warm_start``
        TODO
        See also :term:`partial_fit`.

Attributes
==========

See concept :term:`attribute`.

.. glossary::

    ``classes_``
        TODO

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
