.. _combining_estimators:

===============================================
Pipeline and FeatureUnion: combining estimators
===============================================

.. _pipeline:

Pipeline: chaining estimators
=============================

.. currentmodule:: sklearn.pipeline

:class:`Pipeline` can be used to chain multiple estimators
into one. This is useful as there is often a fixed sequence
of steps in processing the data, for example feature selection, normalization
and classification. :class:`Pipeline` serves two purposes here:

    **Convenience**: You only have to call ``fit`` and ``predict`` once on your
    data to fit a whole sequence of estimators.

    **Joint parameter selection**: You can :ref:`grid search <grid_search>`
    over parameters of all estimators in the pipeline at once.

All estimators in a pipeline, except the last one, must be transformers
(i.e. must have a ``transform`` method).
The last estimator may be any type (transformer, classifier, etc.).


Usage
-----

The :class:`Pipeline` is built using a list of ``(key, value)`` pairs, where
the ``key`` is a string containing the name you want to give this step and ``value``
is an estimator object::

    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.svm import SVC
    >>> from sklearn.decomposition import PCA
    >>> estimators = [('reduce_dim', PCA()), ('clf', SVC())]
    >>> pipe = Pipeline(estimators)
    >>> pipe # doctest: +NORMALIZE_WHITESPACE
    Pipeline(steps=[('reduce_dim', PCA(copy=True, iterated_power='auto',
    n_components=None, random_state=None, svd_solver='auto', tol=0.0,
    whiten=False)), ('clf', SVC(C=1.0, cache_size=200, class_weight=None,
    coef0=0.0, decision_function_shape=None, degree=3, gamma='auto',
    kernel='rbf', max_iter=-1, probability=False, random_state=None,
    shrinking=True, tol=0.001, verbose=False))])

The utility function :func:`make_pipeline` is a shorthand
for constructing pipelines;
it takes a variable number of estimators and returns a pipeline,
filling in the names automatically::

    >>> from sklearn.pipeline import make_pipeline
    >>> from sklearn.naive_bayes import MultinomialNB
    >>> from sklearn.preprocessing import Binarizer
    >>> make_pipeline(Binarizer(), MultinomialNB()) # doctest: +NORMALIZE_WHITESPACE
    Pipeline(steps=[('binarizer', Binarizer(copy=True, threshold=0.0)),
                    ('multinomialnb', MultinomialNB(alpha=1.0,
                                                    class_prior=None,
                                                    fit_prior=True))])

The estimators of a pipeline are stored as a list in the ``steps`` attribute::

    >>> pipe.steps[0]
    ('reduce_dim', PCA(copy=True, iterated_power='auto', n_components=None, random_state=None,
      svd_solver='auto', tol=0.0, whiten=False))

and as a ``dict`` in ``named_steps``::

    >>> pipe.named_steps['reduce_dim']
    PCA(copy=True, iterated_power='auto', n_components=None, random_state=None,
      svd_solver='auto', tol=0.0, whiten=False)

Parameters of the estimators in the pipeline can be accessed using the
``<estimator>__<parameter>`` syntax::

    >>> pipe.set_params(clf__C=10) # doctest: +NORMALIZE_WHITESPACE
    Pipeline(steps=[('reduce_dim', PCA(copy=True, iterated_power='auto',
        n_components=None, random_state=None, svd_solver='auto', tol=0.0,
        whiten=False)), ('clf', SVC(C=10, cache_size=200, class_weight=None,
        coef0=0.0, decision_function_shape=None, degree=3, gamma='auto',
        kernel='rbf', max_iter=-1, probability=False, random_state=None,
        shrinking=True, tol=0.001, verbose=False))])


This is particularly important for doing grid searches::

    >>> from sklearn.model_selection import GridSearchCV
    >>> params = dict(reduce_dim__n_components=[2, 5, 10],
    ...               clf__C=[0.1, 10, 100])
    >>> grid_search = GridSearchCV(pipe, param_grid=params)

Individual steps may also be replaced as parameters, and non-final steps may be
ignored by setting them to ``None``::

    >>> from sklearn.linear_model import LogisticRegression
    >>> params = dict(reduce_dim=[None, PCA(5), PCA(10)],
    ...               clf=[SVC(), LogisticRegression()],
    ...               clf__C=[0.1, 10, 100])
    >>> grid_search = GridSearchCV(pipe, param_grid=params)

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_feature_selection_feature_selection_pipeline.py`
 * :ref:`sphx_glr_auto_examples_model_selection_grid_search_text_feature_extraction.py`
 * :ref:`sphx_glr_auto_examples_plot_digits_pipe.py`
 * :ref:`sphx_glr_auto_examples_plot_kernel_approximation.py`
 * :ref:`sphx_glr_auto_examples_svm_plot_svm_anova.py`

.. topic:: See also:

 * :ref:`grid_search`


Notes
-----

Calling ``fit`` on the pipeline is the same as calling ``fit`` on
each estimator in turn, ``transform`` the input and pass it on to the next step.
The pipeline has all the methods that the last estimator in the pipeline has,
i.e. if the last estimator is a classifier, the :class:`Pipeline` can be used
as a classifier. If the last estimator is a transformer, again, so is the
pipeline.


.. _feature_union:

FeatureUnion: composite feature spaces
======================================

.. currentmodule:: sklearn.pipeline

:class:`FeatureUnion` combines several transformer objects into a new
transformer that combines their output. A :class:`FeatureUnion` takes
a list of transformer objects. During fitting, each of these
is fit to the data independently. For transforming data, the
transformers are applied in parallel, and the sample vectors they output
are concatenated end-to-end into larger vectors.

:class:`FeatureUnion` serves the same purposes as :class:`Pipeline` -
convenience and joint parameter estimation and validation.

:class:`FeatureUnion` and :class:`Pipeline` can be combined to
create complex models.

(A :class:`FeatureUnion` has no way of checking whether two transformers
might produce identical features. It only produces a union when the
feature sets are disjoint, and making sure they are the caller's
responsibility.)


Usage
-----

A :class:`FeatureUnion` is built using a list of ``(key, value)`` pairs,
where the ``key`` is the name you want to give to a given transformation
(an arbitrary string; it only serves as an identifier)
and ``value`` is an estimator object::

    >>> from sklearn.pipeline import FeatureUnion
    >>> from sklearn.decomposition import PCA
    >>> from sklearn.decomposition import KernelPCA
    >>> estimators = [('linear_pca', PCA()), ('kernel_pca', KernelPCA())]
    >>> combined = FeatureUnion(estimators)
    >>> combined # doctest: +NORMALIZE_WHITESPACE
    FeatureUnion(n_jobs=1, transformer_list=[('linear_pca', PCA(copy=True,
        iterated_power='auto', n_components=None, random_state=None,
        svd_solver='auto', tol=0.0, whiten=False)), ('kernel_pca',
        KernelPCA(alpha=1.0, coef0=1, copy_X=True, degree=3,
        eigen_solver='auto', fit_inverse_transform=False, gamma=None,
        kernel='linear', kernel_params=None, max_iter=None, n_components=None,
        n_jobs=1, random_state=None, remove_zero_eig=False, tol=0))],
        transformer_weights=None)


Like pipelines, feature unions have a shorthand constructor called
:func:`make_union` that does not require explicit naming of the components.


Like ``Pipeline``, individual steps may be replaced using ``set_params``,
and ignored by setting to ``None``::

    >>> combined.set_params(kernel_pca=None) # doctest: +NORMALIZE_WHITESPACE
    FeatureUnion(n_jobs=1, transformer_list=[('linear_pca', PCA(copy=True,
          iterated_power='auto', n_components=None, random_state=None,
          svd_solver='auto', tol=0.0, whiten=False)), ('kernel_pca', None)],
        transformer_weights=None)

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_feature_stacker.py`
 * :ref:`sphx_glr_auto_examples_hetero_feature_union.py`
