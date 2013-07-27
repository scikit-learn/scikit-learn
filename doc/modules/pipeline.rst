.. _pipeline:

==============================
Pipeline: chaining estimators
==============================

.. currentmodule:: sklearn.pipeline

:class:`Pipeline` can be used to chain multiple estimators
into one. This is useful as there is often a fixed sequence
of steps in processing the data, for example feature selection, normalization
and classification. :class:`Pipeline` serves two purposes here:

    **Convenience**: You only have to call ``fit`` and ``predict`` once on your
    data to fit a whole sequence of estimators.

    **Joint parameter selection**: You can :ref:`grid search <grid_search>`
    over parameters of all estimators in the pipeline at once.

For estimators to be usable within a pipeline, all except the last one need to have
a ``transform`` function. Otherwise, the dataset can not be passed through this
estimator.


Usage
=====

The :class:`Pipeline` is build using a list of ``(key, value)`` pairs, where
the ``key`` a string containing the name you want to give this step and ``value``
is an estimator object::

    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.svm import SVC
    >>> from sklearn.decomposition import PCA
    >>> estimators = [('reduce_dim', PCA()), ('svm', SVC())]
    >>> clf = Pipeline(estimators)   
    >>> clf # doctest: +NORMALIZE_WHITESPACE
    Pipeline(steps=[('reduce_dim', PCA(copy=True, n_components=None,
        whiten=False)), ('svm', SVC(C=1.0, cache_size=200, class_weight=None,
        coef0=0.0, degree=3, gamma=0.0, kernel='rbf', max_iter=-1,
        probability=False, random_state=None, shrinking=True, tol=0.001,
        verbose=False))])

The estimators of the pipeline are stored as a list in the ``steps`` attribute::

    >>> clf.steps[0]
    ('reduce_dim', PCA(copy=True, n_components=None, whiten=False))

and as a ``dict`` in ``named_steps``::

    >>> clf.named_steps['reduce_dim']
    PCA(copy=True, n_components=None, whiten=False)

Parameters of the estimators in the pipeline can be accessed using the
``<estimator>__<parameter>`` syntax::

    >>> clf.set_params(svm__C=10) # doctest: +NORMALIZE_WHITESPACE
    Pipeline(steps=[('reduce_dim', PCA(copy=True, n_components=None,
        whiten=False)), ('svm', SVC(C=10, cache_size=200, class_weight=None,
        coef0=0.0, degree=3, gamma=0.0, kernel='rbf', max_iter=-1,
        probability=False, random_state=None, shrinking=True, tol=0.001,
        verbose=False))])

This is particularly important for doing grid searches::

    >>> from sklearn.grid_search import GridSearchCV
    >>> params = dict(reduce_dim__n_components=[2, 5, 10],
    ...               svm__C=[0.1, 10, 100])
    >>> grid_search = GridSearchCV(clf, param_grid=params)


.. topic:: Examples:

 * :ref:`example_feature_selection_pipeline.py`
 * :ref:`example_grid_search_text_feature_extraction.py`
 * :ref:`example_plot_digits_pipe.py`
 * :ref:`example_plot_kernel_approximation.py`
 * :ref:`example_svm_plot_svm_anova.py`


Notes
=====

Calling ``fit`` on the pipeline is the same as calling ``fit`` on
each estimator in turn, ``transform`` the input and pass it on to the next step.
The pipeline has all the methods that the last estimator in the pipeline has,
i.e. if the last estimator is a classifier, the :class:`Pipeline` can be used
as a classifier. If the last estimator is a transformer, again, so is the
pipeline.


.. _feature_union:

==========================================
FeatureUnion: Combining feature extractors
==========================================

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
feature sets are disjoint, and making sure they are is the caller's
responsibility.)


Usage
=====

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
        n_components=None, whiten=False)), ('kernel_pca', KernelPCA(alpha=1.0,
        coef0=1, degree=3, eigen_solver='auto', fit_inverse_transform=False,
        gamma=None, kernel='linear', kernel_params=None, max_iter=None,
        n_components=None, remove_zero_eig=False, tol=0))],
        transformer_weights=None)


                                                                       
.. topic:: Examples:

 * :ref:`example_feature_stacker.py`
