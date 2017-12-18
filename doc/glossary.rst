.. currentmodule:: sklearn

.. _glossary:

=========================================
Glossary of Common Terms and API Elements
=========================================

This glossary hopes to definitively represent the tacit and explicit
conventions applied in Scikit-learn and its API, while providing a reference
for users and contributors. It aims to describe the concepts and either detail
their corresponding API or link to other relevant parts of the documentation
which do so.

We begin by listing general concepts (and any that didn't fit elsewhere), but
more specific sets of related terms are listed below:
:ref:`glossary_estimator_types`, :ref:`glossary_target_types`,
:ref:`glossary_methods`, :ref:`glossary_parameters`,
:ref:`glossary_attributes`, :ref:`glossary_sample_props`.

General Concepts
================

.. glossary::

    API
        Refers to both the *specific* interfaces for estimators implemented in
        Scikit-learn and the *generalized* conventions across types of
        estimators as described in this glossary and :ref:`overviewed in the
        contributor documentation <api_overview>`.

        The specific interfaces that constitute Scikit-learn's public API are
        largely documented in :ref:`api_ref`. However we less formally consider
        anything as public API if none of the identifiers required to access it
        begins with ``_``.  We generally try to maintain :term:`backwards
        compatibility` for all objects in the public API.

        Private API, including functions, modules and methods beginning ``_``
        are not assured to be stable.

    array-like
        The most common data format for *input* to Scikit-learn estimators and
        functions, array-like is any type object for which
        :func:`numpy.asarray` will produce an array of appropriate shape
        (usually 1 or 2-dimensional) of appropriate dtype (usually numeric).

        This includes:

        * a numpy array
        * a list of numbers
        * a list of length-k lists of numbers for some fixed length k
        * a :class:`pandas.DataFrame` with all columns numeric
        * a numeric :class:`pandas.Series`

        It excludes:

        * a :term:`sparse matrix`
        * an iterator
        * a generator

        Note that *output* from scikit-learn estimators and functions (e.g.
        predictions) should generally be arrays or sparse matrices, or lists
        thereof (as in multi-output :class:`tree.DecisionTreeClassifier`'s
        ``predict_proba``). An estimator where ``predict()`` returns a list or
        a `pandas.Series` is not valid.

    attribute
    attributes
        We mostly use attribute to refer to how model information is stored on
        an estimator during fitting.  Any public attribute stored on an
        estimator instance is required to begin with an alphabetic character
        and end in a single underscore if it is set in :term:`fit` or
        :term:`partial_fit`.  These are what is documented under an estimator's
        *Attributes* documentation.  The information stored in attributes is
        usually either: sufficient statistics used for prediction or
        transformation; or diagnostic data, such as
        :term:`feature_importances_`.
        Common attributes are listed :ref:`below <glossary_attributes>`.

        Further private attributes used in prediction/transformation/etc. may
        also be set when fitting.  These begin with a single underscore and are
        not assured to be stable for public access.

        A public attribute on an estimator instance that does not end in an
        underscore should be the stored, unmodified value of an ``__init__``
        :term:`parameter` of the same name.  Because of this equivalence, these
        are documented under an estimator's *Parameters* documentation.

    backwards compatibility
        We generally try to maintain backwards compatibility (i.e. interfaces
        and behaviors may be extended but not changed or removed) from release
        to release but this comes with some exceptions:

        Public API only
            The behaviour of objects accessed through private identifiers
            (those beginning ``_``) may be changed arbitrarily between
            versions.
        As documented
            We will generally assume that the users have adhered to the
            documented parameter types and ranges. If the documentation asks
            for a list and the user gives a tuple, we do not assure consistent
            behavior from version to version.
        Deprecation
            Behaviors may change following a :term:`deprecation` period
            (usually two releases long).  Warnings are issued using Python's
            :mod:`warnings` module.
        Keyword arguments
            We may sometimes assume that all optional parameters (other than X
            and y to :term:`fit` and similar methods) are passed as keyword
            arguments only and may be positionally reordered.
        Bug fixes and enhancements
            Bug fixes and -- less often -- enhancements may change the behavior
            of estimators, including the predictions of an estimator trained on
            the same data and :term:`random_state`.  When this happens, we
            attempt to note it clearly in the changelog.
        Serialization
            We make no assurances that pickling an estimator in one version
            will allow it to be unpickled to an equivalent model in the
            subsequent version.  (For estimators in the sklearn package, we
            issue a warning when this unpickling is attempted, even if it may
            happen to work.)  See :ref:`persistence_limitations`.
        :func:`utils.estimator_checks.check_estimator`
            We provide limited backwards compatibility assurances for the
            estimator checks: we may add extra requirements on estimators
            tested with this function, usually when these were informally
            assumed but not formally tested.

        Despite this informal contract with our users, the software is provided
        as is, as stated in the licence.  When a release inadvertently
        introduces changes that are not backwards compatible, these are known
        as software regressions.

    categorical feature
        A categorical or nominal :term:`feature` is one that has a finite set
        of discrete values across the population of data.  These are commonly
        represented as columns of integers or strings. Strings will be rejected
        by most scikit-learn estimators, and integers will be treated as
        ordinal or count-valued. For the use with most estimators, categorical
        variables should be one-hot encoded.
        :class:`~sklearn.preprocessing.CategoricalEncoder` helps encoding
        string-valued categorical features.  Some estimators may handle
        categorical features better when one-hot encoded.  See also
        :ref:`preprocessing_categorical_features` and the
        `http://contrib.scikit-learn.org/categorical-encoding
        <category_encoders>`_ package for tools related to encoding categorical
        features.

    clone
        To copy an :term:`estimator instance` and create a new one with
        identical :term:`parameters`, but without any fitted
        :term:`attributes`, using :func:`~sklearn.base.clone`.

        When ``fit`` is called, a :term:`meta-estimator` usually clones
        a wrapped estimator instance before fitting the cloned instance.
        (Exceptions, for legacy reasons, include
        :class:`~pipeline.Pipeline` and
        :class:`~pipeline.FeatureUnion`.)

    common tests
        This refers to the tests run on almost every estimator class in
        Scikit-learn to check they comply with basic API conventions.  They are
        available for external use through
        :func:`utils.estimator_checks.check_estimator`.

        Note: Some exceptions to the common testing regime are currently
        hard-coded into the library, but we hope to replace this by marking
        exceptional behaviours on the estimator using semantic :term:`estimator
        tags`.

    deprecation
        We use deprecation to slowly violate our :term:`backwards
        compatibility` assurances, usually to to:

        * change the the default value of a parameter; or
        * remove a parameter, attribute, method, class, etc.

        We will ordinarily issue a warning when a deprecated element is used,
        although there may be limitations to this.  For instance, we will raise
        a warning when someone sets a parameter that has been deprecated, but
        may not when they access that parameter's attribute on the estimator
        instance.

        See the :ref:`Contributors' Guide <contributing_deprecation>`.

    docstring
        The embedded documentation for a module, class, function, etc., usually
        in code as a string at the beginning of the object's definition, and
        accessible as the object's ``__doc__`` attribute.

        We try to adhere to `PEP257
        <http://www.python.org/dev/peps/pep-0257/>`_, and follow `NumpyDoc
        conventions <numpydoc.readthedocs.io/en/latest/format.html>`_.

    double underscore
    double underscore notation
        When specifying parameter names for nested estimators, ``__`` may be
        used to separate between parent and child in some contexts. The most
        common use is when setting parameters through a meta-estimator with
        :term:`set_params` and hence in specifying a search grid in
        :ref:`parameter search <grid_search>`. See :term:`parameter`.
        It is also used in :meth:`pipeline.Pipeline.fit` for passing
        :term:`sample properties` to the ``fit`` methods of estimators in
        the pipeline.

    dtype
    data type
        TODO. Mention casting.

    duck typing
        We try to apply `duck typing
        <https://en.wikipedia.org/wiki/Duck_typing>`_ to determine how to
        handle some input values (e.g. checking whether a given estimator is
        a classifier).  That is, we avoid using ``isinstance`` where possible,
        and rely on the presence or absence of attributes to determine an
        object's behaviour.  Some nuance is required when following this
        approach:

        * For some estimators, an attribute may only be available once it is
          :term:`fitted`.  For instance, we cannot a priori determine if
          :term:`predict_proba` is available in a grid search where the grid
          includes alternating between a probabilistic and a non-probabilistic
          predictor in the final step of the pipeline.  In the following, we
          can only determine if ``clf`` is probabilistic after fitting it on
          some data::

              >>> from sklearn.model_selection import GridSearchCV
              >>> from sklearn.linear_model import SGDClassifier
              >>> clf = GridSearchCV(SGDClassifier(),
              ...                    param_grid={'loss': ['log', 'hinge']})

          This means that we can only check for duck-typed attributes after
          fitting, and that we must be careful to make :term:`meta-estimators`
          only present attributes according to the state of the underlying
          estimator after fitting.

        * Checking if an attribute is present (using ``hasattr``) is in general
          just as expensive as getting the attribute (``getattr`` or dot
          notation).  In some cases, getting the attribute may indeed be
          expensive (e.g. for some implementations of
          :term:`feature_importances_`, which may suggest this is an API design
          flaw).  So code which does ``hasattr`` followed by ``getattr`` should
          be avoided; ``getattr`` within a try-except block is preferred.

    estimator instance
        We sometimes use this terminology to distinguish an :term:`estimator`
        class from a constructed instance. For example, in the following,
        ``cls`` is an estimator class, while ``est1`` and ``est2`` are
        instances::

            cls = RandomForestClassifier
            est1 = cls()
            est2 = RandomForestClassifier()

    examples
        We try to give examples of basic usage for most functions and
        classes in the API:

        * as doctests in their docstrings (i.e. within the ``sklearn/`` library
          code itself).
        * as examples in the :ref:`example gallery <general_examples>`
          rendered (using `sphinx-gallery
          <https://sphinx-gallery.readthedocs.io/>`_) from scripts in the
          ``examples/`` directory, exemplifying key features or parameters
          of the estimator/function.
        * sometimes in the :ref:`User Guide <user_guide>` (built from ``doc/``)
          alongside a technical description of the estimator.

    evaluation metric
    evaluation metrics
        TODO

    estimator tags
        A proposed feature (e.g. :issue:`8022`) by which the capabilities of an
        estimator are described through a set of semantic tags.  This would
        enable some runtime behaviors based on estimator inspection, but it
        also allows each estimator to be tested for appropriate invariances
        while being excepted from other :term:`common tests`.

        Some aspects of estimator tags are currently determined through
        the :term:`duck typing` of methods like ``predict_proba`` and through
        some special attributes on estimator objects:

        ``_estimator_type``
            This string-valued attribute identifies an estimator as being a
            classifier, regressor, etc. It is set by mixins such as
            :class:`base.ClassifierMixin`, but needs to be more explicitly
            adopted on a :term:`meta-estimator`.  Its value should usually be
            checked by way of a helper such as :func:`base.is_classifier`.

        ``_pairwise``
            This boolean attribute indicates whether the data (``X``) passed to
            :func:`fit` and similar methods consists of pairwise measures over
            samples rather than a feature representation for each sample.  It
            is usually ``True`` where an estimator has a ``metric`` or
            ``affinity`` or ``kernel`` parameter with value 'precomputed'.
            Its primary purpose is that when a :term:`meta-estimator`
            extracts a sub-sample of data intended for a pairwise estimator,
            the data needs to be indexed on both axes, while other data is
            indexed only on the first axis.

    feature
    features
        In the abstract, a feature is a function (in its mathematical sense)
        mapping a sampled object to a numeric or categorical quantity.
        "Feature" is also commonly used to refer to these quantities, being the
        individual elements of a vector representing a sample. In a data
        matrix, features are represented as columns: each column contains the
        result of applying a feature function to a set of samples.

        Elsewhere features are known as attributes, predictors, regressors, or
        independent variables.

        Nearly all estimators in scikit-learn assume that features are numeric,
        finite and not missing, even when they have semantically distinct
        domains and distributions (categorical, ordinal, count-valued,
        real-valued, interval). See also :term:`categorical feature` and
        :term:`missing values`.

        ``n_features`` indicates the number of features in a dataset.

    fitting
        Calling :term:`fit` (or :term:`fit_transform`, :term:`fit_predict`,
        etc.) on an estimator.

    fitted
        The state of an estimator after :term:`fitting`.

    function
        We provide ad hoc function interfaces for many algorithms, while
        :term:`estimator` classes provide a more consistent interface.

        In particular, Scikit-learn may provide a function interface that fits
        a model to some data and returns the learnt model parameters, as in
        :func:`linear_model.enet_path`.  For transductive models, this also
        returns the embedding or cluster labels, as in
        :func:`manifold.spectral_embedding` or :func:`cluster.dbscan`.  Many
        preprocessing transformers also provide a function interface, akin to
        calling :term:`fit_transform`, as in
        :func:`preprocessing.maxabs_scale`.  Users should be careful to avoid
        :term:`data leakage` when making use of these
        ``fit_transform``-equivalent functions.

        We do not have a strict policy about when to or when not to provide
        function forms of estimators, but maintainers should consider
        consistency with existing interfaces, and whether providing a function
        would lead users astray from best practices (as regards data leakage,
        etc.)

    hyperparameter
    hyper-parameter
        See :term:`parameter`.

    induction
    inductive
        Inductive (contrasted with :term:`transductive`) machine learning
        builds a model of some data that can then be applied to new instances.
        Most estimators in Scikit-learn are inductive, having :term:`predict`
        and/or :term:`transform` methods.

    joblib
        A Python library (http://joblib.readthedocs.io) used in Scikit-learn to
        facilite simple parallelism and caching.  Joblib is oriented towards
        efficiently working with numpy arrays, such as through use of
        :term:`memory mapping`.

    leakage
    data leakage
        A problem in cross validation where generalization performance can be
        over-estimated since knowledge of the test data was inadvertently
        included in training a model.  This is a risk, for instance, when
        applying a :term:`transformer` to the entirety of a dataset rather
        than each training portion in a cross validation split.

        We aim to provide interfaces (such as :mod:`pipeline` and
        :mod:`model_selection`) that shield the user from data leakage.

    memmapping
    memory map
    memory mapping
        TODO

    missing values
        TODO

    ``n_features``
        The number of :term:`features`.

    ``n_outputs``
        The number of :term:`outputs` in the :term:`target`.

    ``n_samples``
        The number of :term:`samples`.

    ``n_targets``
        Synonym for :term:`n_outputs`.

    narrative docs
    narrative documentation
        An alias for :ref:`User Guide <user_guide>`, i.e. documentation written
        in ``doc/modules/``.

    np
        A shorthand for Numpy due to the conventional import statement::

            import numpy as np

    out-of-core
        TODO

    outputs
        TODO?

        Others may call these "responses" or "tasks" or "targets"

    parameter
    parameters
    param
    params
        We mostly use *parameter* to refer to the aspects of an estimator that
        can be specified in its construction. For example, ``max_depth`` and
        ``random_state`` are parameters of :class:`RandomForestClassifier`.
        Parameters to an estimator's constructor are stored unmodified as
        attributes on the estimator instance, and conventionally start with an
        alphabetic character and end with an alphanumeric character.  Each
        estimator's constructor parameters are described in the estimator's
        docstring.

        We do not use parameters in the statistical sense, where parameters are
        values that specify a model and can be estimated from data. What we
        call parameters might be what statisticians call hyperparameters to the
        model: aspects for configuring model structure that are often not
        directly learnt from data.  However, our parameters are also used to
        prescribe modeling operations that do not affect the learnt model, such
        as :term:`n_jobs` for controlling parallelism.

        When talking about the parameters of a :term:`meta-estimator`, we may
        also be including the parameters of the estimators wrapped by the
        meta-estimator.  Ordinarily, these nested parameters are denoted by
        using a :term:`double underscore` (``__``) to separate between the
        estimator-as-parameter and its parameter.  Thus ``clf =
        BaggingClassifier(base_estimator=DecisionTreeClassifier(max_depth=3))``
        has a deep parameter ``base_estimator__max_depth`` with value ``3``,
        which is accessible with ``clf.base_estimator.max_depth`` or
        ``clf.get_params()['base_estimator__max_depth']``.

        The list of parameters and their current values can be retrieved from
        an :term:`estimator instance` using its :term:`get_params` method.

        Between construction and fitting, parameters may be modified using
        :term:`set_params`.  To enable this, parameters are not ordinarily
        validated or altered when the estimator is constructed, or when each
        parameter is set. Parameter validation is performed when :term:`fit` is
        called.

        Common parameters are listed :ref:`below <glossary_parameters>`.

    pairwise metric
        TODO

        See precomputed.

    pd
        A shorthand for `Pandas <http://pandas.pydata.org>`_ due to the
        conventional import statement::

            import pandas as pd

    precomputed
        TODO

    rectangular
        Data that can be represented as a matrix with :term:`samples` on the
        first axis and a fixed, finite set of :term:`features` on the second
        is called rectangular.

        This term excludes samples with non-vectorial structure, such as text,
        an image of arbitrary size, a time series of arbitrary length, a set of
        vectors, etc. The purpose of a :term:`vectorizer` is to produce
        rectangular forms of such data.

    sample
    samples
        We usually use this term as a noun to indicate a single feature vector.
        Elsewhere a sample is called an instance, data point, or observation.
        ``n_samples`` indicates the number of samples in a dataset, being the
        number of rows in a data array ``X``.

    sample property
    sample properties
        A sample property is data for each sample (e.g. an array of length
        n_samples) passed to an estimator method or a similar function,
        alongside but distinct from the :term:`features` (``X``) and
        :term:`target` (``y``). The most prominent example is
        :term:`sample_weight`; see others at :ref:`glossary_sample_props`.

        As of version 0.19 we do not have a consistent approach to handling
        sample properties and their routing in :term:`meta-estimators`, though
        a ``fit_params`` parameter is often used.

    scikit-learn-contrib
        A venue for publishing Scikit-learn-compatible libraries that are
        broadly authorized by the core developers and the contrib community,
        but not maintained by the core developer team.
        See http://scikit-learn-contrib.github.io.

    sparse matrix
        A representation of two-dimensional numeric data that is more memory
        efficient the corresponding dense numpy array where almost all elements
        are zero. We use the :mod:`scipy.sparse` framework, which provides
        several underlying sparse data representations, or *formats*.
        Some formats are more efficient than others for particular tasks, and
        when a particular format provides especial benefit, we try to document
        this fact in Scikit-learn parameter descriptions.

        Some sparse matrix formats (notably CSR, CSC, COO and LIL) distinguish
        between *implicit* and *explicit* zeros. Explicit zeros are stored
        (i.e. they consume memory in a ``data`` array) in the data structure,
        while implicit zeros correspond to every element not otherwise defined
        in explicit storage.

        Two semantics for sparse matrices are used in Scikit-learn:

        matrix semantics
            The sparse matrix is interpreted as an array with implicit and
            explicit zeros being interpreted as the number 0.  This is the
            interpretation most often adopted, e.g. when sparse matrices
            are used for feature matrices or multilabel indicator matrices.
        graph semantics
            As with :mod:`scipy.sparse.csgraph`, explicit zeros are
            interpreted as the number 0, but implicit zeros indicate a masked
            or absent value, such as the absence of an edge between two
            vertices of a graph, where an explicit value indicates an edge's
            weight. This interpretation is adopted to represent connectivity
            in clustering, in representations of nearest neighborhoods
            (e.g. :func:`neighbors.kneighbors_graph`), and for precomputed
            distance representation where only distances in the neighborhood
            of each point are required.

        When working with sparse matrices, we assume that it is sparse for a
        good reason, and avoid writing code that densifies a user-provided
        sparse matrix, instead maintaining sparsity or raising an error if not
        possible (i.e. if an estimator does not / cannot support sparse
        matrices).

    target
    targets
        TODO

        Dependent variable or outcome variable.

    transduction
    transductive
        A transductive (contrasted with :term:`inductive`) machine learning
        method is designed to model a specific dataset, but not to apply that
        model to unseen data.  Examples include :class:`manifold.TSNE`,
        :class:`cluster.AgglomerativeClustering` and
        :class:`neighbors.LocalOutlierFactor`.

    unlabeled
    unlabeled data
        TODO

.. _glossary_estimator_types:

Class APIs and Estimator Types
==============================

.. glossary::

    classifier
        TODO
        Mention that within scikit-learn, all support multi-class
        classification, defaulting to OvR.
        Mention :func:`~base.is_classifier`.

    clusterer
        TODO

    estimator
        TODO

        The core functionality of some estimators may also be available as a
        :term:`function`.

    feature extractor
    feature extractors
        A :term:`transformer` which takes input where each sample is not
        represented as an :term:`array-like` object of fixed length, and
        produces an `array-like` object of :term:`features` for each sample
        (and thus a 2-dimensional array-like for a set of samples).  In other
        words, it (lossily) maps a non-rectangular data representation into
        :term:`rectangular` data. Feature extractors should define
        :term:`get_feature_names`.

    meta-estimator
    meta-estimators
        TODO

        Mention duck typing. Mention that duck typing of methods only works
        after fitting (see
        :func:`utils.metaestimators.if_delegate_has_method`). Mention lenient
        validation. ?Mention sample props.  Mention clone.

    outlier detector
        TODO

    predictor
        An estimator which provides :term:`predict`.  This encompasses
        :term:`classifier`, :term:`regressor`, :term:`outlier detector` and
        sometimes :term:`clusterer` (although the latter only provide
        :term:`fit_predict` when they are not :term:`inductive`).

    regressor
        TODO
        Mention :func:`~base.is_regressor`.

    transformer
        An estimator supporting :term:`fit` and :term:`transform` and/or
        :term:`fit_transform`.  A purely :term:`transductive` transformer,
        such as :class:`manifold.TSNE`, may not implement ``transform``.

    vectorizer
        See :term:`feature extractor`.

    ``Xt``
        Shorthand for transformed ``X``.

There are further APIs specifically related to a small family of estimators,
such as:

.. glossary:

    cross validation splitter
    CV splitter
        A non-estimator family of classes used to split a dataset into a
        sequence of train and test portions (see :ref:`cross_validation`),
        by providing :term:`split` and :term:`get_n_splits` methods.
        Note that unlike estimators, these do not have :term:`fit` methods
        and do not provide :term:`set_params` or :term:`get_params`.
        Parameter validation may be performed in ``__init__``.

    scorer
        A non-estimator callable object which evaluates an estimator on given
        test data, returning a number. See :ref:`scoring_parameter`; see also
        :term:`evaluation metric`.

Further examples:

* :class:`neighbors.DistanceMetric`
* :class:`gaussian_process.kernels.Kernel`
* ``tree.Criterion``

.. _glossary_target_types:

Target Types
============

.. glossary::

    binary
        A classification problem consisting of two classes.  A binary target
        may represented as for a :term:`multiclass` problem but with only two
        labels.  A binary decision function is represented as a 1d array.

        Semantically, one class is often considered the "positive" class.
        Unless otherwise specified (e.g. using :term:`pos_label` in
        :term:`evaluation metrics`), we consider the class label with the
        greater value (numerically or lexicographically) as the positive class:
        of labels [0, 1], 1 is the positive class; of [1, 2], 2 is the positive
        class; of ['no', 'yes'], 'yes' is the positive class; of ['no', 'YES'],
        'no' is the positive class.  This affects the output of
        :term:`decision_function`, for instance.

        Note that a dataset sampled from a multiclass ``y`` or a continuous
        ``y`` may appear to be binary.

        :func:`~utils.multiclass.type_of_target` will return 'binary' for
        binary input, or a similar array with only a single class present.

    continuous
        TODO

    multiclass
        A classification problem consisting of more than two classes.  A
        multiclass target may be represented as a 1-dimensional array of
        strings or integers.  A 2d column vector of integers (i.e. a
        single output in :term:`multioutput` terms) is also accepted.

        We do not officially support other orderable, hashable objects as class
        labels, even if estimators may happen to work when given classification
        targets of such type.

        For semi-supervised classification, :term:`unlabeled` samples should
        have the special label -1 in ``y``.

        Within sckit-learn, all estimators supporting binary classification
        also support multiclass classification, using One-vs-Rest by default.

        A :class:`preprocessing.LabelEncoder` helps to canonicalize multiclass
        targets as integers.

        :func:`~utils.multiclass.type_of_target` will return 'multiclass' for
        multiclass input. The user may also want to handle 'binary' input
        identically to 'multiclass'.

    multilabel
        A :term:`multioutput` target where each output is :term:`binary`.  This
        may be represented as a 2d (dense) array or sparse matrix of integers,
        such that each column is a separate binary target, where positive
        labels are indicated with 1 and negative labels are usually -1 or 0.
        Sparse multilabel targets are not supported everywhere that dense
        multilabel targets are supported.

        Semantically, a multilabel target can be thought of as a set of labels
        for each sample.  While not used internally,
        :class:`preprocessing.MultiLabelBinarizer` is provided as a utility to
        convert from a list of sets representation to a 2d array or sparse
        matrix. One-hot encoding a multiclass target with
        :class:`preprocessing.LabelBinarizer` turns it into a multilabel
        problem.

        :func:`~utils.multiclass.type_of_target` will return
        'multilabel-indicator' for multilabel input, whether sparse or dense.

    multioutput continuous
        TODO

    multioutput multiclass
        TODO

        A classification problem...

        :mod:`multioutput` provides estimators which estimate multi-output
        problems using multiple single-output estimators.  This may not fully
        account for dependencies among the different outputs, which methods
        natively handling the multioutput case (e.g. decision trees, nearest
        neighbors, neural networks) may do better.

.. _glossary_methods:

Methods
=======

.. glossary::

    ``decision_function``
        In a fitted :term:`classifier` or :term:`outlier detector`, predicts a
        "soft" score for each sample in relation to each class, rather than the
        "hard" categorical prediction produced by :term:`predict`.  Its input
        is usually only some observed data, :term:`X`.

        Output conventions:

        binary classification
            A 1-dimensional array, where values strictly greater than zero
            indicate the positive class (i.e. the last class in
            :term:`classes_`).
        multiclass classification
            A 2-dimensional array, where the row-wise arg-maximum is the
            predicted class.  Columns are ordered according to
            :term:`classes_`.
        multilabel classification
            Scikit-learn is inconsistent in its representation of multilabel
            decision functions.  Some estimators represent it like multioutput
            multiclass, i.e. a list of 2d arrays, each with two columns. Others
            represent it with a single 2d array, whose columns correspond to
            the individual binary classification decisions. The latter
            representation is ambiguously identical to the multiclass
            classification format, though its semantics differ: it should be
            interpreted, like in the binary case, by thresholding at 0.

            TODO
            ..
            see https://gist.github.com/jnothman/4807b1b0266613c20ba4d1f88d0f8cf5
        multioutput classification
            A list of 2d arrays, corresponding to each multiclass decision
            function.
        outlier detection
            A 1-dimensional array, where a value greater than or equal to zero
            (TODO: check equality case) indicates an inlier.

    ``fit``
        The ``fit`` method is provided on every estimator. It usually takes some
        :term:`samples` ``X``, :term:`targets` ``y`` if the model is supervised,
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

        :ref:`glossary_target_types` describes possible formats for ``y``.

    ``fit_predict``
        Used especially for :term:`transductive` estimators, this fits the
        model and returns the predictions (similar to :term:`predict`) on the
        training data. In clusterers, these predictions are also stored in the
        :term:`labels_` attribute, and the output of ``.fit_predict(X)`` is
        usually equivalent to ``.fit(X).predict(X)``. The parameters to
        ``fit_predict`` are the same as those to ``fit``.

    ``fit_transform``
        A method on :term:`transformers` which fits the estimator and returns
        the transformed training data. It takes parameters as in :term:`fit`
        and its output should have the same shape as calling ``.fit(X,
        ...).transform(X)``. There are nonetheless rare cases where
        ``.fit_transform(X, ...)`` and ``.fit(X, ...).transform(X)`` do not
        return the same value, wherein training data needs to be handled
        differently (due to model blending in stacked ensembles, for instance;
        such cases should be clearly documented).
        :term:`Transductive` transformers may also provide ``fit_transform``
        but not :term:`transform`.

        One reason to implement ``fit_transform`` is that performing ``fit``
        and ``transform`` separately would be less efficient than together.
        :class:`base.TransformerMixin` provides a default implementation,
        providing a consistent interface across transformers where
        ``fit_transform`` is or is not specialised.

        In :term:`inductive` learning -- where the goal is to learn a
        generalised model that can be applied to new data -- users should be
        careful not to apply ``fit_transform`` to the entirety of a dataset
        (i.e. training and test data together) before further modelling, as
        this results in :term:`data leakage`.

    ``get_feature_names``
        Primarily for :term:`feature extractors`, but also used for other
        transformers to provide string names for each column in the output of
        the estimator's :term:`transform` method.  It outputs a list of
        strings, and may take a list of strings as input, corresponding
        to the names of input columns from which output column names can
        be generated.  By default input features are named x0, x1, ....

    ``get_n_splits``
        On a :term:`CV splitter` (not an estimator), returns the number of
        elements one would get if iterating through the return value of
        :term:`split` given the same parameters.  Takes the same parameters as
        split.

    ``get_params``
        Gets all :term:`parameters`, and their values, that can be set using
        :term:`set_params`.  A parameter ``deep`` can be used, when set to
        False to only return those parameters not including ``__``, i.e.  not
        due to indirection via contained estimators.
        
        Most estimators adopt the definition from :class:`base.BaseEstimator`,
        which simply adopts the parameters defined for ``__init__``.
        :class:`pipeline.Pipeline`, among others, reimplements ``get_params``
        to declare the estimators named in its ``steps`` parameters as
        themselves being parameters.

    ``partial_fit``
        Facilitates fitting an estimator in an online fashion.  Unlike ``fit``,
        repeatedly calling ``partial_fit`` does not clear the model, but
        updates it with respect to the data provided. The portion of data
        provided to ``partial_fit`` may be called a mini-batch.
        Each mini-batch must be of consistent shape, etc.

        ``partial_fit`` may also be used for :term:`out-of-core` learning,
        although limited to the case where learning can be performed online,
        i.e. the model is usable after each ``partial_fit`` and there is no
        separate processing needed to finalize the model.
        :class:`cluster.Birch` introduces the convention that calling
        ``partial_fit(X)`` will produce a model that is not finalized, but the
        model can be finalized by calling ``partial_fit()`` i.e. without
        passing a further mini-batch.

        Generally, estimator parameters should not be modified between calls
        to ``partial_fit``, although ``partial_fit`` should validate them
        as well as the new mini-batch of data.  In contrast, ``warm_start``
        is used to repeatedly fit the same estimator with the same data
        but varying parameters.

        Like ``fit``, ``partial_fit`` should return the estimator object.

        To clear the model, a new estimator should be constructed, for instance
        with :func:`base.clone`.

    ``predict``
        Makes a prediction for each sample, usually only taking :term:`X` as
        input (but see under regressor output conventions below). In a
        :term:`classifier` or :term:`regressor`, this prediction is in the same
        target space used in fitting (e.g. one of {'red', 'amber', 'green'} if
        the ``y`` in fitting consisted of these strings).  Despite this, even
        when ``y`` passed to :term:`fit` is a list or other array-like, the
        output of ``predict`` should always be an array or sparse matrix. In a
        :term:`clusterer` or :term:`outlier detector` the prediction is an
        integer.

        Output conventions:

        classifier
            An array of shape ``(n_samples,)`` ``(n_samples, n_outputs)``.
            :term:`Multilabel` data may be represented as a sparse matrix if
            a sparse matrix was used in fitting. Each element should be one
            of the values in the classifier's :term:`classes_` attribute.

        clusterer
            An array of shape ``(n_samples,) where each value is from 0 to
            ``n_clusters - 1`` if the corresponding sample is clustered,
            and -1 if the sample is not clustered, as in
            :func:`cluster.dbscan`.

        outlier detector
            An array of shape ``(n_samples,)`` where each value is -1 for an
            outlier and 1 otherwise.

        regressor
            A numeric array of shape ``(n_samples,)``, usually float64.
            Some regressors have extra options in their ``predict`` method,
            allowing them to return standard deviation (``return_std=True``)
            or covariance (``return_cov=True``) relative to the predicted
            value.  In this case, the return value is a tuple of arrays
            corresponding to (prediction mean, std, cov) as required.

    ``predict_log_proba``
        The natural logarithm of the output of :term:`predict_proba`, provided
        to facilitate numerical stability.

    ``predict_proba``
        A method in classifiers that are able to return probability estimates
        for each class.  Its input is usually only some observed data,
        :term:`X`.

        Output conventions are like those for :term:`decision_function` except
        in the :term:`binary` classification case, where one column is output
        for each class (while ``decision_function`` outputs a 1d array). For
        binary and multiclass predictions, each row should add to 1.

        Like other methods, ``predict_proba`` should only be present when the
        estimator can make probabilistic predictions (see :term:`duck typing`).
        This means that the presence of the method may depend on estimator
        parameters (e.g. in :class:`linear_model.SGDClassifier`) or training
        data (e.g. in :class:`model_selection.GridSearchCV`) and may only
        appear after fitting.

    ``score``
        A method on an estimator, usually a :term:`predictor`, which evaluates
        its predictions on a given dataset, and returns a single numerical
        score.  A greater return value should indicate better predictions;
        accuracy is used for classifiers and R^2 for regressors by default.

    ``score_samples``
        TODO

    ``set_params``
        Available in any estimator, takes keyword arguments corresponding to
        keys in :term:`get_params`.  Each is provided a new value to assign
        such that calling ``get_params`` after ``set_params`` will reflect the
        changed :term:`parameters`.  Most estimators use the implementation in
        :class:`base.BaseEstimator`, which handles nested parameters and
        otherwise sets the parameter as an attribute on the estimator.
        The method is overridden in :class:`pipeline.Pipeline` and related
        estimators.

    ``split``
        On a :term:`CV splitter` (not an estimator), this method accepts
        parameters (:term:`X`, :term:`y`, :term:`groups`), where all may be
        optional, and returns an iterator over ``(train_idx, test_idx)``
        pairs.  Each of {train,test}_idx is a 1d integer array, with values
        from 0 from ``X.shape[0] - 1`` of any length, such that no values
        appear in both some ``train_idx`` and its corresponding ``test_idx``.

    ``transform``
        In a :term:`transformer`, transforms the input, usually only :term:`X`,
        into some transformed space (conventionally notated as :term:`Xt`).
        Output is an array or sparse matrix of length :term:`n_samples` and
        with number of columns fixed after :term:`fitting`.

.. _glossary_parameters:

Parameters
==========

These common parameter names, specifically used in estimator construction
(see concept :term:`parameter`), sometimes also appear as parameters of
functions or non-estimator constructors.

.. glossary::

    ``class_weight``
        Used to specify sample weights when fitting classifiers as a function
        of the :term:`target` class.  Where :term:`sample_weight` is also
        supported and given, it is multiplied by the ``class_weight``
        contribution. Similarly, where ``class_weight`` is used in a
        :term:`multioutput` (including :term:`multilabel`) tasks, the weights
        are multiplied across outputs (i.e. columns of ``y``).

        By default all samples have equal weight such that classes are
        effectively weighted by their their prevalence in the training data.
        This could be achieved explicitly with ``class_weight={label1: 1,
        label2: 1, ...}`` for all class labels.

        More generally, ``class_weight`` is specified as a dict mapping class
        labels to weights (``{class_label: weight}``), such that each sample
        of the named class is given that weight.

        ``class_weight='balanced'`` can be used to give all classes
        equal weight by giving each sample a weight inversely related
        to its class's prevalence in the training data:
        ``n_samples / (n_classes * np.bincount(y))``.

        For multioutput classification, a list of dicts is used to specify
        weights for each output. For example, for four-class multilabel
        classification weights should be ``[{0: 1, 1: 1}, {0: 1, 1: 5}, {0: 1,
        1: 1}, {0: 1, 1: 1}]`` instead of ``[{1:1}, {2:5}, {3:1}, {4:1}]``.

        The ``class_weight`` parameter is validated and interpreted with
        :func:`utils.compute_class_weight`.

    ``cv``
        Determines a cross validation splitting strategy, as used in
        cross-validation based routines. ``cv`` is also available in estimators
        such as :class:`multioutput.ClassifierChain` or
        :class:`calibration.CalibratedClassifierCV` which use the predictions
        of one estimator as training data for another, to not overfit the
        training supervision.

        Possible inputs for ``cv`` are usually:

        - An integer, specifying the number of folds in K-fold cross
          validation. K-fold will be stratified over classes if the estimator
          is a classifier (determined by :func:`base.is_classifier`) and the
          :term:`targets` may represent a binary or multiclass (but not
          multioutput) classification problem (determined by
          :func:`utils.multiclass.type_of_target`).
        - A :term:`cross validation splitter` instance. Refer to the
          :ref:`User Guide <cross_validation>` for splitters available
          within Scikit-learn.
        - An iterable yielding train/test splits.

        With some exceptions (especially where not using cross validation at
        all is an option), the default is 3-fold.

        ``cv`` values are validated and interpreted with :func:`utils.check_cv`.

    ``kernel``
        TODO

    ``max_iter``
        For estimators involving iterative optimization, this determines the
        maximum number of iterations to be performed in :term:`fit`.  If
        ``max_iter`` iterations are run without convergence, a
        :class:`exceptions.ConvergenceWarning` should be raised.  Note that the
        interpretation of "a single iteration" is inconsistent across
        estimators: some, but not all, use it to mean a single epoch (i.e. a
        pass over every sample in the data).

        FIXME perhaps we should have some common tests about the relationship
        between ConvergenceWarning and max_iter.

    ``memory``
        Some estimators make use of :class:`joblib.Memory` to
        store partial solutions during fitting. Thus when ``fit`` is called
        again, those partial solutions have been memoized and can be reused.

        A ``memory`` parameter can be specified as a string with a path to a
        directory, or a :class:`joblib.Memory` instance (or an object with a
        similar interface, i.e. a ``cache`` method) can be used.

        ``memory`` values are validated and interpreted with
        :func:`utils.validation.check_memory`.

    ``metric``
        As a parameter, this is the scheme for determining the distance between
        two data points.  See :func:`metrics.pairwise_distances`.  In practice,
        for some algorithms, an improper distance metric (one that does not
        obey the triangle inequality, such as Cosine Distance) may be used.

        XXX: hierarchical clustering uses ``affinity`` with this meaning.

        We also use *metric* to refer to :term:`evaluation metrics`, but avoid
        using this sense as a parameter name.

    ``n_components``
        TODO

    ``n_jobs``
        This is used to specify how many concurrent processes/threads should be
        used for parallelized routines.  Scikit-learn uses one processor for
        its processing by default, although it also makes use of NumPy, which
        may be configured to use a threaded numerical processor library (like
        MKL; see :ref:`FAQ <faq_mkl_threading>`).

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
          thing.  *It can be highly detrimental to performance to run multiple
          copies of some estimators or functions in parallel.*

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
                number of different distinct random seeds. Popular integer
                random seeds are 0 and `42
                <http://en.wikipedia.org/wiki/Answer_to_the_Ultimate_Question_of_Life%2C_the_Universe%2C_and_Everything>`_.

            A :class:`numpy.random.RandomState` instance
                Use the provided random state, only affecting other users
                of the same random state instance. Calling fit multiple times
                will reuse the same instance, and will produce different
                results.

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
        estimator's :term:`score` method is used.  See :ref:`scoring_parameter`
        in the User Guide.

        Where multiple metrics can be evaluated, ``scoring`` may be given
        either as a list of unique strings or a dict with names as keys and
        callables as values. Note that this does *not* specify which score
        function is to be maximised, and another parameter such as ``refit``
        may be used for this purpose.

        The ``scoring`` parameter is validated and interpreted using
        :func:`metrics.check_scoring`.

    ``verbose``
        Logging is not handled very consistently in Scikit-learn at present,
        but when it is provided as an option, the ``verbose`` parameter is
        usually available to choose no logging (set to False). Any True value
        should enable some logging, but larger integers (e.g. above 10) may be
        needed for full verbosity.  Verbose logs are usually printed to
        Standard Output.
        Estimators should not produce any output on Standard Output with the
        default ``verbose`` setting.

    ``warm_start``

        When fitting an estimator repeatedly on the same dataset, but for
        multiple parameter values (such as to find the value maximizing
        performance as in :ref:`grid search <grid_search>`), it may be possible
        to reuse aspects of the model learnt from the previous parameter value,
        saving time.  When ``warm_start`` is true, the existing :term:`fitted`
        model :term:`attributes` an are used to initialise the new model
        in a subsequent call to :term:`fit`.

        Note that this is only applicable for some models and some
        parameters, and even some orders of parameter values. For example,
        ``warm_start`` may be used when building random forests to add more
        trees to the forest (increasing ``n_estimators``) but not to reduce
        their number.

        :term:`partial_fit` also retains the model between calls, but differs:
        with ``warm_start`` the parameters change and the data is constant
        across calls to ``fit``; with ``partial_fit``, the mini-batch of data
        changes and model parameters stay fixed.

.. _glossary_attributes:

Attributes
==========

See concept :term:`attribute`.

.. glossary::

    ``classes_``
        A list of class labels known to the :term:`classifier`, mapping each
        label to a numerical index used in the model representation our output.
        For instance, the array output from :term:`predict_proba` has columns
        aligned with ``classes_``. For :term:`multi-output` classifiers,
        ``classes_`` should be a list of lists, with one class listing for
        each output.  For each output, the classes should be sorted
        (numerically, or lexicographically for strings).

        ``classes_`` and the mapping to indices is often managed with
        :class:`preprocessing.LabelEncoder`.

    ``components_``
        An affine transformation matrix of shape ``(n_components, n_features)``
        used in many linear :term:`transformers` where :term:`n_components` is
        the number of output features and :term:`n_features`` is the number of
        input features.

        See also :term:`components_` which is a similar attribute for linear
        predictors.

    ``coef_``
        The weight/coefficient matrix of a generalised linear model
        :term:`predictor`, of shape ``(n_features,)`` for binary classification
        and single-output regression, ``(n_classes, n_features)`` for
        multiclass classification and ``(n_targets, n_features)`` for
        multi-output regression. Note this does not include the intercept
        (or bias) term, which is stored in ``intercept_``.

        When available, ``feature_importances_`` is not usually provided as
        well, but can be calculated as the  norm of each feature's entry in
        ``coef_``.

        See also :term:`components_` which is a similar attribute for linear
        transformers.

    ``embedding_``
        An embedding of the training data in :ref:`manifold learning
        <manifold>` estimators, with shape ``(n_samples, n_components)``,
        identical to the output of :term:`fit_transform`.  See also
        :term:`labels_`.

    ``n_iter_``
        The number of iterations actually performed when fitting an iterative
        estimator that may stop upon convergence. See also :term:`max_iter`.

    ``feature_importances_``
        A vector of shape ``(n_features,)`` available in some
        :term:`predictors` to provide a relative measure of the importance of
        each feature in the predictions of the model.

    ``labels_``
        A vector containing a cluster label for each sample of the training
        data in :term:`clusterers`, identical to the output of
        :term:`fit_predict`.  See also :term:`embedding_`.

.. _glossary_sample_props:

Sample properties
=================

See concept :term:`sample property`.

.. glossary::

    ``groups``
        Used in cross validation routines to identify samples which are
        correlated.  Each value is an identifier such that, in a supporting
        :term:`CV splitter`, samples from some ``groups`` value may not
        appear in both a training set and its corresponding test set.
        See :ref:`group_cv`.

    ``sample_weight``
        A relative weight for each sample.  Intuitively, if all weights are
        integers, a weighted model or score should be equivalent to that
        calculated when repeating the sample the number of times specified in
        the weight.  Weights may be specified as floats, so that sample weights
        are usually equivalent up to a constant positive scaling factor.

        FIXME  Is this interpretation always the case in practice? We have no
        common tests.

        Some estimators, such as decision trees, support negative weights.
        FIXME: This feature or its absence may not be tested or documented in
        many estimators.

        This is not entirely the case where other parameters of the model
        consider the number of samples in a region, as with ``min_samples`` in
        :class:`cluster.DBSCAN`.  In this case, a count of samples becomes
        to a sum of their weights.

        In classification, sample weights can also be specified as a function
        of class with the :term:`class_weight` estimator :term:`parameter`.

    ``X``
        Denotes data that is observed at training and prediction time, used as
        independent variables in learning.  The notation is uppercase to denote
        that it is ordinarily a matrix (see :term:`rectangular`).

    ``Xt``
        Shorthand for "transformed :term:`X`".

    ``y``
    ``Y``
        Denotes data that may be observed at training time as the dependent
        variable in learning, but which is unavailable at prediction time, and
        is usually the :term:`target` of prediction.  The notation may be
        uppercase to denote that it is a matrix, representing
        :term:`multi-output` targets, for instance; but usually we use ``y``
        and sometimes do so even when multiple outputs are assumed.

