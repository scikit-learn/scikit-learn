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
        coef0=0.0, degree=3, gamma=0.0, kernel='rbf', probability=False,
        shrinking=True, tol=0.001, verbose=False))])

The estimators of the pipeline are stored as a list in the ``steps`` attribute::

    >>> clf.steps[0]
    ('reduce_dim', PCA(copy=True, n_components=None, whiten=False))

and as a ``dict`` in ``named_steps``::

    >>> clf.named_steps['reduce_dim']
    PCA(copy=True, n_components=None, whiten=False)

Parameters of the estimators in the pipeline can be accessed using the
``<estimator>__<parameter>`` syntax::

    >>> clf.set_params(svm__C=10) # NORMALIZE_WHITESPACE
    Pipeline(steps=[('reduce_dim', PCA(copy=True, n_components=None, whiten=False)), ('svm', SVC(C=10, cache_size=200, class_weight=None, coef0=0.0, degree=3, gamma=0.0,
      kernel='rbf', probability=False, shrinking=True, tol=0.001,
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
The pipeline has all the methods that the last estimator in the pipline has,
i.e. if the last estimator is a classifier, the :class:`Pipeline` can be used
as a classifier. If the last estimator is a transformer, again, so is the
pipeline.
